from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
from astropy.table import Table
import os, sys, copy, pickle, numpy as np

basedir='/home/ckilpatrick/scripts/python/progenitors/'
params = {
    'SHEET':'1paDfeYsJyv9X_XL26gV9Pk2xF9VjG70T4Wly5lRPuxs',
    'token':basedir+'token.pickle',
    'credentials':basedir+'credentials.json',
    'target': basedir,
    'yse': {'user': 'ckilpatrick', 'password': 'Vfg190OW@K9E*g4$Bpmw'},
    'cols': ['Name',   'YSEPZ',   'TNS', 'OSC', 'RA',  'Dec',
        'Classification', 'Host',    'NED',
        'Discovery Date','HST','Spectrum','Distance','Distance Method',
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

        return(alldata)

def convert_table_to_lists(table):

    output = [params['cols']]

    for row in table:
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
            table_data = convert_table_to_lists(all_tables[key])
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
                          'endIndex': 199,},
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

            # Mask if adding mask
            if mask and 'mask' in all_tables[key].meta.keys():
                reqs = []
                for idx in np.where(all_tables[key].meta['mask'])[0]:
                    req={
                      'updateDimensionProperties': {
                        'range': {
                          'sheetId': sheetid,
                          'dimension': 'ROWS',
                          'startIndex': int(idx)+1,
                          'endIndex': int(idx)+2,
                        },
                        'properties': {
                          'hiddenByUser': True,
                        },
                        'fields': 'hiddenByUser',
                    }}
                    reqs.append(req)

                if len(reqs)>0:
                    body = {'requests': reqs}
                    dum = sheet.batchUpdate(spreadsheetId=spreadsheetId,
                      body=body).execute()

def load_metadata_from_file(all_data, directory):
    for key in all_data.keys():
        savekey = key.replace(' ','_').replace('-','_').replace('/','_')
        file = directory+savekey+'.pkl'
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

def save_metadata_to_file(all_data, directory):
    for key in all_data.keys():
        savekey = key.replace(' ','_').replace('-','_').replace('/','_')
        with open(directory+savekey+'.pkl','wb') as f:
            pickle.dump(all_data[key].meta, f)

