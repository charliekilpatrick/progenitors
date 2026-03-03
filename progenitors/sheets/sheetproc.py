"""Google Sheets API integration for progenitor data download/upload.

Settings: sheet column names and OAuth scopes in progenitors/settings/sheets.py
"""
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
from astropy.table import Table
import os
import sys
import copy
import pickle
import numpy as np

from ..settings.sheets import SHEET_COLUMNS, SCOPES

# Repo root: from progenitors/sheets/sheetproc.py go up to repo
_here = os.path.abspath(os.path.dirname(__file__))
basedir = os.path.dirname(os.path.dirname(_here))
_credentials = os.environ.get("PROGENITORS_CREDENTIALS_PATH") or os.path.join(basedir, "credentials.json")

params = {
    "SHEET": os.environ.get("PROGENITORS_SHEET", ""),
    "token": os.path.join(basedir, "token.pickle"),
    "credentials": _credentials,
    "target": basedir,
    "yse": {
        "user": os.environ.get("YSE_USER", ""),
        "password": os.environ.get("YSE_PASSWORD", ""),
    },
    "cols": list(SHEET_COLUMNS),
}


def transpose(data):
    return list(map(list, zip(*data)))


def initiate_gsheet(token_file, credentials=None):
    creds = None
    if os.path.exists(token_file):
        with open(token_file, 'rb') as token:
            creds = pickle.load(token)
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(credentials, SCOPES)
            creds = flow.run_local_server(port=0)
        with open(token_file, 'wb') as token:
            pickle.dump(creds, token)
    service = build('sheets', 'v4', credentials=creds)
    return service.spreadsheets()


def download_progenitor_data(spreadsheetId):
    sheet = initiate_gsheet(params['token'], credentials=params['credentials'])
    response = sheet.get(spreadsheetId=spreadsheetId).execute()
    if not response or 'sheets' not in response:
        return {}
    alldata = {}
    for subsheet in response['sheets']:
        name = subsheet['properties']['title']
        data = sheet.values().get(spreadsheetId=spreadsheetId, range=name).execute()
        if 'values' in data:
            header = data['values'][0]
            bodydata = data['values'][1:]
            bodydata = [b for b in bodydata if len(b[0].strip()) > 0]
            tabdata = transpose(bodydata)
            if len(tabdata) == 0:
                tabdata = [['']] * len(header)
            elif len(tabdata) < len(header):
                for _ in range(len(header) - len(tabdata)):
                    tabdata.append([' ' * 100] * len(tabdata[0]))
            table = Table(tabdata, names=header)
            table.meta['name'] = name
            remove_rows = [i for i, row in enumerate(table) if len(row['Name'].strip()) == 0]
            table.remove_rows(remove_rows)
            alldata[name] = table
        else:
            init_data = [['X' * 100] for _ in params['cols']]
            table = Table(init_data, names=params['cols'])
            alldata[name] = table
    return alldata


def convert_table_to_lists(table, mask=None):
    output = [list(params['cols'])]
    for i, row in enumerate(table):
        if mask is not None and mask[i]:
            continue
        add_row = []
        for k in output[0]:
            if k in row.colnames:
                add_row.append('\'' + str(row[k]))
            elif k.lower() in row.colnames:
                add_row.append('\'' + str(row[k.lower()]))
            else:
                add_row.append('')
        output.append(add_row)
    return output


def upload_progenitor_data(spreadsheetId, all_tables, mask=False):
    sheet = initiate_gsheet(params['token'])
    response = sheet.get(spreadsheetId=spreadsheetId).execute()
    if not response or 'sheets' not in response:
        return
    for key in all_tables.keys():
        mask_vals = all_tables[key].meta.get('mask') if mask else None
        table_data = convert_table_to_lists(all_tables[key], mask=mask_vals)
        sheetid_list = [s['properties']['sheetId'] for s in response['sheets']
                       if s['properties']['title'] == key]
        if len(sheetid_list) == 0:
            body = {'requests': [{'addSheet': {'properties': {'title': key}}}]}
            sheet.batchUpdate(spreadsheetId=spreadsheetId, body=body).execute()
            newresp = sheet.get(spreadsheetId=spreadsheetId).execute()
            sheetid = [s['properties']['sheetId'] for s in newresp['sheets']
                     if s['properties']['title'] == key][0]
        else:
            sheetid = sheetid_list[0]
        body = {'requests': [{'updateDimensionProperties': {
            'range': {'sheetId': sheetid, 'dimension': 'ROWS', 'startIndex': 0, 'endIndex': 9999},
            'properties': {'hiddenByUser': False},
            'fields': 'hiddenByUser',
        }}]}
        sheet.batchUpdate(spreadsheetId=spreadsheetId, body=body).execute()
        body = {'requests': [{'updateCells': {'range': {'sheetId': sheetid}, 'fields': 'userEnteredValue'}}]}
        sheet.batchUpdate(spreadsheetId=spreadsheetId, body=body).execute()
        body = {'majorDimension': 'ROWS', 'values': table_data}
        sheet.values().update(spreadsheetId=spreadsheetId, range=key, valueInputOption='USER_ENTERED', body=body).execute()


def load_metadata_from_file(all_data, directory):
    for key in all_data.keys():
        savekey = key.replace(' ', '_').replace('-', '_').replace('/', '_')
        file = os.path.join(directory, savekey + '.pkl')
        print('Loading metadata from:', file)
        if os.path.exists(file):
            try:
                with open(file, 'rb') as f:
                    all_data[key].meta = pickle.load(f)
            except (EOFError, MemoryError) as e:
                print(f'Error opening {file}, skipping...', e)
    return all_data


def save_metadata_to_file(all_data, directory, make_backup=False):
    for key in all_data.keys():
        savekey = key.replace(' ', '_').replace('-', '_').replace('/', '_')
        os.makedirs(directory, exist_ok=True)
        fulloutname = os.path.join(directory, savekey + '.pkl')
        backupname = os.path.join(directory, 'backup', savekey + '.backup.pkl')
        print(f'Saving {key} data to {fulloutname}')
        with open(fulloutname, 'wb') as f:
            pickle.dump(all_data[key].meta, f)
        if make_backup:
            os.makedirs(os.path.join(directory, 'backup'), exist_ok=True)
            with open(backupname, 'wb') as f:
                pickle.dump(all_data[key].meta, f)
