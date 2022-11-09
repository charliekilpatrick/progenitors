import util
import sheetproc
import sys
from astropy.table import unique

def main(redo=False, download_yse=True, update_classification=False):

    # Build out progenitor data from google sheet and add metadata from files
    sndata = sheetproc.download_progenitor_data(util.params['SHEET'])
    sndata = sheetproc.load_metadata_from_file(sndata, util.params['metadata'])

    # Add new targets from YSEPZ if there are any
    if download_yse: sndata = util.add_yse_targets(sndata)

    # Iterate through data by spectroscopic type - Ia, IIb, IIn, II, Ib/c, Other
    all_keys = list(sndata.keys())
    for type_key in all_keys:

        for row in sndata[type_key]:
            all_spectra = util.get_spectra(row, sndata[type_key])

        # Update type_data with new information from OSC, TNS, YSEPZ, NED
        for dattype in ['ned','distance','tns','osc','yse','ads','hst']:
            print(redo)
            update = util.add_metadata(sndata[type_key], dattype, redo=redo)
            if vars(update)!=vars(sndata[type_key]):
                sndata[type_key]=copy.copy(update)
                sheetproc.save_metadata_to_file({type_key: sndata[type_key]},
                    util.params['metadata'])

        # Finally add column data to the table
        for coltype in ['OSC','TNS','YSEPZ','Host','Discovery Date',
            'Classification','Distance','Redshift','NED','HST',
            'Distance Method','Ref. (Distance)',
            'Ref. (Classification)','Ref. (Discovery)']:
            print('Running',coltype,'for',type_key)
            sndata[type_key] = util.add_data(sndata[type_key], coltype)

        # If there are any entries in the table, make sure they're unique
        if len(sndata[type_key])>1:
            sndata[type_key] = unique(sndata[type_key], keys='Name')

        # Sort final table by Discovery Date
        if 'Discovery Date' in sndata[type_key].keys():
            sndata[type_key].sort('Discovery Date')

        sndata[type_key].meta['mask']=util.get_classification_mask(
            sndata[type_key])

        # Upload the table data to google sheets and save metadata to a file
        sheetproc.upload_progenitor_data(util.params['SHEET'], sndata,
            mask=True)
        sheetproc.save_metadata_to_file(sndata, util.params['metadata'])

    # Rearrange tables based on classifications
    if update_classification:
        new_sndata = {}
        for key in all_keys:
            print('Updating classifications',key)
            new_sndata[key] = util.gather_type(sndata, key)

            print('Objects in sheet',key,len(new_sndata[key]))

            # Sort final table by Discovery Date
            if 'Discovery Date' in new_sndata[key].keys():
                new_sndata[key].sort('Discovery Date')

            new_sndata[key].meta['mask']=util.get_classification_mask(
                new_sndata[key])

        sheetproc.upload_progenitor_data(util.params['SHEET'], new_sndata,
            mask=True)
        sheetproc.save_metadata_to_file(new_sndata, util.params['metadata'])

main(redo=False, update_classification=True)
