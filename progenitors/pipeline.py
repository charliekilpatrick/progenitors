"""
Main pipeline: build progenitor data from Google Sheets and metadata sources.

Downloads sheet data, fetches metadata (NED, TNS, YSE, ADS, OSC, MAST, etc.),
optionally updates classification and sends email alerts. Settings (default
metadata list, trim keys, YSE query) are in progenitors.settings.pipeline.
"""
import sys
import os
import time
import copy
import pickle
from astropy.table import unique, Table

from . import util
from .sheets import sheetproc
from . import options
from .settings.pipeline import (
    TRIM_METADATA_KEEP_KEYS,
    DEFAULT_METADATA,
    DEFAULT_YSE_SQL_QUERY,
)


def do_trim_metadata(sndata, keep_keys=TRIM_METADATA_KEEP_KEYS):
    """
    Remove metadata entries whose object name is not in the table.

    Parameters
    ----------
    sndata : dict of Table
        Per-sheet tables with .meta dicts.
    keep_keys : set or list, optional
        Metadata keys to always keep (e.g. 'all_sndata').

    Returns
    -------
    dict
        sndata (modified in place).
    """
    for key in sndata.keys():
        init_keys = list(sndata[key].meta.keys())
        for obj in init_keys:
            if obj in keep_keys:
                continue
            if obj not in sndata[key]['Name'].data:
                print(f'Deleting {obj} from {key} metadata')
                del sndata[key].meta[obj]
    return sndata


def main(
    redo=False,
    download_yse=True,
    update_classification=False,
    alert=True,
    yse_sql_query=DEFAULT_YSE_SQL_QUERY,
    always_update=False,
    trim_metadata=False,
    redo_obj=None,
    metadata=DEFAULT_METADATA,
):
    """
    Run the full progenitor pipeline: Sheets, metadata, optional classification and alerts.

    Parameters
    ----------
    redo : bool, optional
        Re-fetch metadata even if already present.
    download_yse : bool, optional
        Add YSE targets to sheet data.
    update_classification : bool, optional
        Reclassify objects across sheets and upload.
    alert : bool, optional
        Send email for new candidates when updating classification.
    yse_sql_query : str, optional
        SQL for YSE target query.
    always_update : bool, optional
        Always save metadata after each type.
    trim_metadata : bool, optional
        Call do_trim_metadata before final save.
    redo_obj : list of str, optional
        Object names to redo metadata for only.
    metadata : list of str, optional
        Metadata types to fetch (e.g. 'ned', 'distance', 'tns').
    """
    sndata = sheetproc.download_progenitor_data(util.params['SHEET'])
    sndata = sheetproc.load_metadata_from_file(sndata, util.params['metadata'])

    for key in sndata.keys():
        if 'all_sndata' in sndata[key].meta.keys():
            meta = copy.copy(sndata[key].meta)
            curr_table = Table(sndata[key], meta={})
            all_sn_data = Table(sndata[key].meta['all_sndata'], meta={})
            sndata[key] = all_sn_data
            sndata[key].meta = meta
            sndata[key].meta['curr_table'] = curr_table

    if download_yse:
        sndata = util.add_yse_targets(sndata, yse_sql_query=yse_sql_query)

    for key in sndata.keys():
        sndata[key].meta['mask'] = util.get_classification_mask(sndata[key])

    all_keys = list(sndata.keys())
    for type_key in all_keys:
        meta_size = sys.getsizeof(pickle.dumps(sndata[type_key].meta)) / (1024.0 ** 2)
        meta_size = float('%.3f' % meta_size)
        print(f'Current metadata file size {meta_size} MB')

        for dattype in metadata:
            update = util.add_metadata(sndata[type_key], dattype, redo=redo, redo_obj=redo_obj)
            if vars(update) != vars(sndata[type_key]) or always_update:
                print(f'Updating {type_key} metadata in pickle file')
                sndata[type_key] = copy.copy(update)
                sheetproc.save_metadata_to_file({type_key: sndata[type_key]}, util.params['metadata'])

        for coltype in sheetproc.params['cols']:
            print('Running', coltype, 'for', type_key)
            sndata[type_key] = util.add_data(sndata[type_key], coltype)

        if len(sndata[type_key]) > 1:
            sndata[type_key] = unique(sndata[type_key], keys='Name')
        sndata[type_key].sort('Discovery Date')
        sndata[type_key].meta['mask'] = util.get_classification_mask(sndata[type_key])

        sheetproc.upload_progenitor_data(util.params['SHEET'], sndata, mask=True)
        new_size = sys.getsizeof(pickle.dumps(sndata[type_key].meta)) / (1024.0 ** 2)
        new_size = float('%.3f' % new_size)
        if new_size > meta_size:
            print(f'New metadata {new_size} MB > {meta_size} MB. Saving...')
            sheetproc.save_metadata_to_file({type_key: sndata[type_key]}, util.params['metadata'])

    if trim_metadata:
        do_trim_metadata(sndata)
    sheetproc.save_metadata_to_file(sndata, util.params['metadata'], make_backup=True)

    if update_classification:
        new_sndata = {}
        for key in all_keys:
            print('Updating classifications', key)
            new_sndata[key] = util.gather_type(sndata, key)
            print('Objects in sheet', key, len(new_sndata[key]))

        for key in all_keys:
            for row in sndata[key]:
                if any(row['Name'] in new_sndata[k]['Name'] for k in new_sndata):
                    continue
                new_sndata['Other'].add_row(row)
                new_sndata['Other'].meta[row['Name']] = sndata[key].meta[row['Name']]

        for key in new_sndata.keys():
            new_sndata[key] = unique(new_sndata[key], keys='Name')
            new_sndata[key].sort('Discovery Date')
            new_sndata[key].meta['mask'] = util.get_classification_mask(new_sndata[key])
            if alert and 'curr_table' in new_sndata[key].meta.keys():
                mask = new_sndata[key].meta['mask']
                curr_table = new_sndata[key].meta['curr_table']
                for ii, row in enumerate(new_sndata[key]):
                    if not mask[ii] and row['Name'] not in curr_table['Name'].data:
                        util.post_alert(row)

        sheetproc.upload_progenitor_data(util.params['SHEET'], new_sndata, mask=True)
        for key in new_sndata.keys():
            copy_table_data = copy.copy(new_sndata[key])
            copy_table_data.meta = None
            new_sndata[key].meta['all_sndata'] = copy_table_data
        sheetproc.save_metadata_to_file(new_sndata, util.params['metadata'], make_backup=True)
