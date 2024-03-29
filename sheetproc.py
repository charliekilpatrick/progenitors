from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
from astropy.table import Table
import os, sys, copy, pickle, numpy as np

basedir=os.path.split(os.path.realpath(__file__))[0]
params = {
    'SHEET':os.environ['PROGENITORS_SHEET'],
    'token':os.path.join(basedir, 'token.pickle'),
    'credentials':os.path.join(basedir, 'credentials.json'),
    'target': basedir,
    'yse': {'user': os.environ['YSE_USER'],
            'password': os.environ['YSE_PASSWORD']},
    'cols': ['Name',   'YSEPZ',   'TNS', 'RA',  'Dec',
        'Classification', 'Host',    'NED',
        'Discovery Date','HST (pre-explosion)','HST (post-explosion)',
        'JWST (pre-explosion)','JWST (post-explosion)',
        'Spectrum','Distance','Distance Method',
        'Ref. (Distance)', 'Redshift',    'Ref. (Discovery)',
        'Ref. (Classification)',   'Post-Explosion']
}

SCOPES = ['https://www.googleapis.com/auth/spreadsheets']

def transpose(data):
    return(list(map(list, zip(*data))))

def initiate_gsheet(token_file, credentials=None):
    creds = None
    if os.path.exists(token_file):
        with open(token_file, 'rb') as token:
            creds = pickle.load(token)
    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(
                credentials, SCOPES)
            creds = flow.run_local_server(port=0)
        # Save the credentials for the next run
        with open(token_file, 'wb') as token:
            pickle.dump(creds, token)

    service = build('sheets', 'v4', credentials=creds)
    sheet = service.spreadsheets()

    return(sheet)

def download_progenitor_data(spreadsheetId):
    sheet = initiate_gsheet(params['token'], credentials=params['credentials'])
    response = sheet.get(spreadsheetId=spreadsheetId).execute()

    if response and 'sheets' in response.keys():
        alldata = {}

        for subsheet in response['sheets']:
            name=subsheet['properties']['title']
            subsheetId=subsheet['properties']['sheetId']

            data = sheet.values().get(spreadsheetId=spreadsheetId,
                range=name).execute()

            if 'values' in data.keys():
                header = data['values'][0]
                bodydata = data['values'][1:]
                bodydata = [b for b in bodydata if len(b[0].strip())>0]
                tabdata = transpose(bodydata)
                if len(tabdata)==0:
                    tabdata = [['']]*len(header)
                elif len(tabdata)<len(header):
                    for i in np.arange(len(header)-len(tabdata)):
                        newcol = [' '*100]*len(tabdata[0])
                        tabdata.append(newcol)

                table = Table(tabdata, names=header)
                table.meta['name'] = name

                remove_rows = []
                for i,row in enumerate(copy.copy(table)):
                    if len(row['Name'].strip())==0:
                        remove_rows.append(i)

                table.remove_rows(remove_rows)

                alldata[name] = table

            else:
                init_data = []
                header = []
                for key in params['cols']:
                    header.append(key)
                    init_data.append(['X'*100])

                newtable = Table(init_data, names=header)

                alldata[name] = newtable

        return(alldata)

def convert_table_to_lists(table, mask=None):

    output = [params['cols']]

    for i,row in enumerate(table):
        if mask is not None and mask[i]:
            continue
        add_row = []
        for k in output[0]:
            if k in row.colnames: add_row.append('\''+str(row[k]))
            elif k.lower() in row.colnames: add_row.append('\''+str(row[k.lower()]))
            else: add_row.append('')
        output.append(add_row)

    return(output)

def upload_progenitor_data(spreadsheetId, all_tables, mask=False):
    sheet = initiate_gsheet(params['token'])
    response = sheet.get(spreadsheetId=spreadsheetId).execute()

    if response and 'sheets' in response.keys():

        for key in all_tables.keys():
            mask_vals = None
            if mask and 'mask' in all_tables[key].meta.keys():
                mask_vals = all_tables[key].meta['mask']

            table_data = convert_table_to_lists(all_tables[key], mask=mask_vals)

            sheetid = [s['properties']['sheetId'] for s in response['sheets']
                if s['properties']['title']==key]

            if len(sheetid)==0:
                body = {
                    'requests': [{
                        'addSheet': {
                            'properties': {
                                'title': key,
                            }
                        }
                    }]
                }
                req = sheet.batchUpdate(spreadsheetId=spreadsheetId,
                    body=body).execute()
                newresp = sheet.get(spreadsheetId=spreadsheetId).execute()
                sheetid = [s['properties']['sheetId'] for s in newresp['sheets']
                    if s['properties']['title']==key][0]
            else:
                sheetid=sheetid[0]

            # Unhide all rows
            body={'requests':[{'updateDimensionProperties': {
                        'range': {
                          'sheetId': sheetid,
                          'dimension': 'ROWS',
                          'startIndex': 0,
                          'endIndex': 9999,},
                          'properties': {'hiddenByUser': False,},
                          'fields': 'hiddenByUser',}}]}
            dum = sheet.batchUpdate(spreadsheetId=spreadsheetId,
                  body=body).execute()

            # Clear all previous values
            body={'requests':[{'updateCells': {
                        'range': {'sheetId': sheetid},
                        'fields': 'userEnteredValue'}}]}
            dum = sheet.batchUpdate(spreadsheetId=spreadsheetId,
                  body=body).execute()

            # Update data to table_data
            body = {'majorDimension': 'ROWS', 'values': table_data}
            result = sheet.values().update(spreadsheetId=spreadsheetId,
                range=key, valueInputOption='USER_ENTERED', body=body).execute()

def load_metadata_from_file(all_data, directory):
    for key in all_data.keys():
        savekey = key.replace(' ','_').replace('-','_').replace('/','_')
        file = os.path.join(directory, savekey+'.pkl')
        print('Loading metadata from:',file)
        if os.path.exists(file):
            try:
                with open(file, 'rb') as f:
                    all_data[key].meta = pickle.load(f)
            except EOFError:
                print('Error opening {0}, skipping...'.format(file))
            except MemoryError:
                print('Error opening {0}, skipping...'.format(file))

    return(all_data)

def save_metadata_to_file(all_data, directory, make_backup=False):
    for key in all_data.keys():
        savekey = key.replace(' ','_').replace('-','_').replace('/','_')
        if not os.path.exists(directory):
            os.makedirs(directory)
        fulloutname = os.path.join(directory, savekey+'.pkl')
        backupname = os.path.join(directory, 'backup', savekey+'.backup.pkl')
        print(f'Saving {key} data to {fulloutname}')
        with open(fulloutname, 'wb') as f:
            pickle.dump(all_data[key].meta, f)
        if make_backup:
            if not os.path.exists(os.path.join(directory, 'backup')):
                os.makedirs(os.path.join(directory, 'backup'))
            with open(backupname, 'wb') as f:
                pickle.dump(all_data[key].meta, f)

