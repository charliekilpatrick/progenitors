import util
import sheetproc
import sys
from astropy.table import unique, Table
import os
import time
import glob
import shutil
import copy
import options
import pickle

def do_trim_metadata(sndata, keep_keys=['all_sndata', 'curr_table','mask',
    'jwst','hst']):

    for key in sndata.keys():
        init_keys = list(sndata[key].meta.keys())
        for obj in init_keys:
            # Don't delete additional data for table versions or mask
            if obj in keep_keys: continue
            if obj not in sndata[key]['Name'].data:
                # Basically, get rid of metadata for an object if it is
                # not currently in the full table for this type.  This can
                # happen, e.g., when an object is reclassified but its
                # metadata remain in the previous table
                print(f'Deleting {obj} from {key} metadata')
                del sndata[key].meta[obj]

    return(sndata)


def main(redo=False, download_yse=True, update_classification=False,
    alert=True, yse_sql_query='160', always_update=False, trim_metadata=False,
    redo_obj=None, metadata=['osc','jwst','ned','distance','tns',
        'yse','ads','hst']):

    # Build out progenitor data from google sheet and add metadata from files
    sndata = sheetproc.download_progenitor_data(util.params['SHEET'])
    sndata = sheetproc.load_metadata_from_file(sndata, util.params['metadata'])

    for key in sndata.keys():
        if 'all_sndata' in sndata[key].meta.keys():
            # Need to save metadata so it is not overwritten
            meta = copy.copy(sndata[key].meta)
            # Now create current data table and recover table with all data
            curr_table = Table(sndata[key], meta={})
            all_sn_data = Table(sndata[key].meta['all_sndata'], meta={})
            # Overwrite the current table with all data so that we don't
            # re-analyze objects that are already in the table
            sndata[key] = all_sn_data
            sndata[key].meta = meta
            sndata[key].meta['curr_table']=curr_table

    # Add new targets from YSEPZ if there are any
    if download_yse:
        sndata = util.add_yse_targets(sndata, yse_sql_query=yse_sql_query)

    # Do initial update of classification masks to reflect the current number
    # and ordering of objects in the table
    for key in sndata.keys():
        sndata[key].meta['mask']=util.get_classification_mask(sndata[key])

    # Iterate through data by spectroscopic type - Ia, IIb, IIn, II, Ib/c, Other
    all_keys = list(sndata.keys())
    for type_key in all_keys:
        # Get size of metadata to see if we need to resave
        meta_size = sys.getsizeof(pickle.dumps(sndata[type_key].meta))
        meta_size = meta_size / (1024.) / (1024.)
        meta_size = float('%.3f'%meta_size)
        print(f'Current metadata file size {meta_size} MB')

        # Update type_data with new information from TNS, YSEPZ, NED
        for dattype in metadata:
            update = util.add_metadata(sndata[type_key], dattype, redo=redo,
                redo_obj=redo_obj)
            if vars(update)!=vars(sndata[type_key]) or always_update:
                print(f'Updating {type_key} metadata in pickle file')
                sndata[type_key]=copy.copy(update)
                # This will only save the current type key
                sheetproc.save_metadata_to_file({type_key: sndata[type_key]},
                    util.params['metadata'])

        # Finally add column data to the table
        # See params['cols'] in sheetproc.py for type and ordering of columns
        for coltype in sheetproc.params['cols']:
            print('Running',coltype,'for',type_key)
            sndata[type_key] = util.add_data(sndata[type_key], coltype)

        # If there are any entries in the table, make sure they're unique
        if len(sndata[type_key])>1:
            sndata[type_key] = unique(sndata[type_key], keys='Name')

        # Sort final table by Discovery Date
        sndata[type_key].sort('Discovery Date')

        # Update classification mask after new metadata
        sndata[type_key].meta['mask']=util.get_classification_mask(
            sndata[type_key])

        # Upload the table data to google sheets and save metadata to a file
        sheetproc.upload_progenitor_data(util.params['SHEET'], sndata,
            mask=True)
        new_size = sys.getsizeof(pickle.dumps(sndata[type_key].meta))
        new_size = new_size / (1024.) / (1024.)
        new_size = float('%.3f'%new_size)
        # Only save if there's been an appreciable addition of data since we're
        # going to save at the end anyway.  This is mainly a check for instances
        # where we added a lot of data to avoid loss due to an error.
        if new_size > meta_size:
            print(f'New metadata {new_size} MB > {meta_size} MB. Saving...')
            sheetproc.save_metadata_to_file({type_key: sndata[type_key]},
                    util.params['metadata'])

    # Trim metadata before saving files if we want to perform this check
    if trim_metadata: do_trim_metadata(sndata)

    # Always do a final save
    sheetproc.save_metadata_to_file(sndata, util.params['metadata'],
        make_backup=True)

    # Rearrange tables based on classifications
    if update_classification:
        new_sndata = {}
        for key in all_keys:
            print('Updating classifications',key)
            new_sndata[key] = util.gather_type(sndata, key)

            print('Objects in sheet',key,len(new_sndata[key]))

        # If there are any objects not added to new_sndata, add to 'Other'
        for key in all_keys:
            for row in sndata[key]:
                if any([row['Name'] in new_sndata[key]['Name']
                    for key in new_sndata.keys()]):
                    continue

                new_sndata['Other'].add_row(row)
                new_sndata['Other'].meta[row['Name']]=sndata[key].meta[row['Name']]

        for key in new_sndata.keys():
            new_sndata[key] = unique(new_sndata[key], keys='Name')
            new_sndata[key].sort('Discovery Date')

            # Need to do this step last so mask matches with ordering of table
            new_sndata[key].meta['mask']=util.get_classification_mask(
                new_sndata[key])

            if alert and 'curr_table' in new_sndata[key].meta.keys():
                mask = new_sndata[key].meta['mask']
                curr_table = new_sndata[key].meta['curr_table']
                for ii,row in enumerate(new_sndata[key]):
                    if (not mask[ii] and
                        row['Name'] not in curr_table['Name'].data):
                        util.post_alert(row)

        sheetproc.upload_progenitor_data(util.params['SHEET'], new_sndata,
            mask=True)

        for key in new_sndata.keys():
            copy_table_data = copy.copy(new_sndata[key])
            copy_table_data.meta = None
            new_sndata[key].meta['all_sndata'] = copy_table_data

        # Always do a final save
        sheetproc.save_metadata_to_file(new_sndata, util.params['metadata'],
            make_backup=True)

if __name__ == '__main__':
    start = time.time()

    args = options.parse_arguments()

    command = ' '.join(sys.argv)

    if args.metadata is not None:
        metadata = args.metadata.split(',')
    else:
        metadata = ['osc','jwst','ned','distance','tns','yse','ads','hst']

    # YSE sql queries - 160=just z<0.02 (but has to have redshift), 378=all HST
    main(redo=args.redo, update_classification=args.update_classification,
          alert=args.alert, yse_sql_query=args.yse_sql_query,
          always_update=args.always_update, trim_metadata=args.trim_metadata,
          redo_obj=args.redo_obj, metadata=metadata)

    total_time = time.time()-start
    message = f'Finished with: {command}\n'
    message += f'It took {total_time} seconds to complete this script.'
    options.message(message)
