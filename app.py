#!/usr/bin/env python
# coding: utf-8

# In[1]:


import base64
import urllib.parse
import gzip, io
from time import time

import pandas as pd
from scipy.spatial import distance_matrix
import umap
from sklearn.preprocessing import StandardScaler

import dash
from dash.dependencies import Input, Output, State
from dash import dcc
from dash import html
from dash import dash_table
import plotly.graph_objs as go
import plotly.express as px


# In[2]:


pythonScriptDownloadLink = 'assets/In-Silico_Screener_VCF_merger.py'
linkStyle = {'width' : '20%', 'margin' : 10, 'display' : 'inline-block', 'color' : "#077be2",'font-family' : 'gisha', 'fontSize' : 15, 'textAlign' : 'center' ,'borderWidth' : '2px','borderColor' : "#000044",'borderStyle' : 'groove','padding': 4,'borderRadius' : '5px'}


# In[3]:


external_stylesheets = ['assets/ISS.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server
app.title = "In-Silico Screener | Ohad Birk's Lab"
#app.config['suppress_callback_exceptions'] = True


# In[4]:


app.layout = html.Div(
    dcc.Loading(
    type = 'dot',
    children = [
    
    html.Div([ #header div
        html.Img(id = 'header_image',
            src = 'assets/logo_birklab.png',
            style = {
                'height' : '20%',
                'width' : '20%',
                'float' : 'left'
            }
        ),
        html.Span(''),
        html.Img(id = 'header_image2',
            src = 'assets/in-silico-screener.png',
            style = {
                'height' : '20%',
                'width' : '20%',
                'float' : 'right',
                'margin' : 4
            }
        ),
    ],
    style={'width' : '100%','textAlign' : 'center'}
    ),
    
    html.Div([ #upload Div
        html.P(
            dcc.Upload(
                id='yyf_upload',
                children=html.Div(
                    [
                    'Drag or ',
                    html.A('select YYF file here')
                    ],
                    id = 'yyf_upload_text',
                    style={
                        'textAlign' : 'center', 'fontSize' : 18, 'padding': 12, 'borderRadius' : '5px', 'width':'40%', 'marginBottom' : 140,
                        'font-family' : 'gisha', 'lineHeight' : '100%', 'borderWidth' : '2px',
                        'borderColor' : '#000044', 'borderStyle' : 'dashed', 'marginLeft' : 'auto', 'marginRight' : 'auto'
                    }
                ),
            style={'height':'60px','marginLeft':'auto','marginRight':'auto','width':'100%'},
            multiple = False
            ),
            style={'width': '100%', 'display': 'inline-block'}
        ),

        html.Br(),

        
        html.P(
            dcc.Upload(
                id='groups_upload',
                children=html.Div([
                'optional: add ',
                html.A('groups file')
                ],
                id = 'grouping_upload_text',
                style={
                    'textAlign' : 'center', 'fontSize' : 18, 'padding': 12, 'borderRadius' : '5px', 'width':'40%',
                    'font-family' : 'gisha', 'lineHeight' : '100%', 'borderWidth' : '2px',
                    'borderColor' : '#000044', 'borderStyle' : 'dashed', 'marginLeft' : 'auto', 'marginRight' : 'auto'
                    }
                ),
                style={'height':'60px','marginLeft':'auto','marginRight':'auto','width':'100%'},
                    multiple = False
            ),
            style={'width': '100%', 'display': 'inline-block'}
        ),
        
        html.Div([
            html.Button(
                'Analyze',
                id = 'Analysis_start_btn',
                n_clicks = 0,
                disabled = True,
                style = {'text-transform' : 'none', 'font-size' : 24, 'margin' : 'auto'},
                ),
            ],
            style = {'width' : '100%', 'display' : 'inline-block','textAlign' : 'center'}
        ),
        html.Div([
            html.P('Awating input files...',
                   id = 'upload_status',
                   style = {'font-size' : 16, 'width' : '100%', 'display' : 'inline-block','textAlign' : 'center'}
                ),
            ],
            style = {'width' : '100%', 'display' : 'inline-block','textAlign' : 'center'}
        ),
        
        html.Div([
            html.P(html.A('Download Python Script',target='_blank', href = pythonScriptDownloadLink, style = {'fontSize' : '120%'}), style = linkStyle),
            html.Br(),
            html.Iframe(
                src = 'assets/Instructions.pdf',
                style={'height': '600px', 'width' : '80%', 'display' : 'inline-block','textAlign' : 'center'}
            ),
            ],
            style = {'width' : '80%', 'display' : 'inline-block','textAlign' : 'center'}
        ),
        
        ], id = 'upload_div',
        ), #end of upload div
        
        
        html.Div([
            dash_table.DataTable(
                id = 'table',
                filter_action = 'native',
                sort_action = 'native',
                fixed_rows={'headers' : True},
                style_cell = {
                    'textAlign': 'left',
                    'overflow': 'hidden',
                    'textOverflow': 'ellipsis',
                    'maxWidth': 0,
                    'font_size': 15
                },
                tooltip_duration=None,
                style_data_conditional=[
                    {
                    'if': {'row_index': 'odd'},
                    'backgroundColor': 'rgb(250, 250, 255)'
                    },
                    {
                    'if': {'row_index': 'even'},
                    'backgroundColor': 'rgb(252, 252, 255)'
                    },
                ],
                style_cell_conditional=[
                    {'if': {'column_id' : 'ClinVar'},
                     'color': 'tomato',
                     'fontWeight': 'bold'
                    },
                ],
                style_data={
                    'whiteSpace': 'normal',
                    'height': 'auto',
                    },
                style_header = {
                    'backgroundColor': 'rgb(220, 220, 255)',
                    'fontWeight': 'bold',
                    'whiteSpace' : 'normal'
                    },
                style_table={
                    'maxHeight': '250px',
                    'maxwidth' : '150%',
                    'border': 'thin lightgrey solid'
                    },
                style_as_list_view = False
            ),

            html.Div([
                html.Div(
                    id = 'information_window_div',
                    style = {'marginTop' : 30, 'marginLeft' : 30, 'width' : '100%', 'display' : 'inline-block', 'textAlign' : 'left'}
                ),
                ], style = {'width' : '100%', 'display' : 'inline-block', 'textAlign' : 'center'}
            ),
            html.Hr(),
            html.Button(
                'Check sample similarities',
                id = 'population_start_btn',
                n_clicks = 0,
                style = {'text-transform' : 'none', 'font-size' : 24, 'margin' : 'auto'},
            ),
            html.Div([
                html.Div(
                    id = 'population_analysis_div',
                    style = {'marginTop' : 10, 'marginLeft' : 0, 'width' : '100%', 'display' : 'inline-block', 'textAlign' : 'left'}
                ),
                ], style = {'width' : '100%', 'display' : 'inline-block', 'textAlign' : 'center'}
            ),
        ],
        style = {'display' : 'none'},
        id = 'results_div',
        ),
        dcc.Store(id = 'referenceGenome'),
        dcc.Store(id = 'yyf_df'),
        dcc.Store(id = 'grouping_df'),
        dcc.Store(id = 'heatmap_df'),
    ]
    )
)


# In[5]:


@app.callback(
    [
    Output('Analysis_start_btn', 'disabled'),
    Output('yyf_upload_text', 'children'),
    Output('yyf_df', 'data'),
    ],
    [
    Input('yyf_upload', 'contents'),
    ],
    [
    State('yyf_upload', 'filename'),
    ]
)
def uploadYYF(yyf, file_name):
    if yyf != None:
        try:
            content_type, content_string = yyf.split(',')
            decoded = base64.b64decode(content_string)
            if file_name.lower().endswith('.gz'):
                decoded = io.BytesIO(decoded)
                gzf = gzip.GzipFile(fileobj=decoded)
                decoded = gzf.read().decode('utf-8-sig', 'ignore')
            else:
                decoded = decoded.decode('utf-8-sig', 'ignore')
            data = [line.strip().split(',') for line in decoded.split("\n") if len(line) > 2]
            if data[0][:3] != ['Alleles','Phenotypes','Gene'] or data[0][3] not in ('Coordinates - hg19', 'Coordinates - hg38') or data[0][4:7] != ['ClinVar ID', 'Severity', 'Review Status']:
                return [False, ['Corrupted YYF file, ', html.A('reupload')], None]
            return [False,['YYF file uploaded'], data]
        except:
            return [True, ['Corrupted YYF file, ', html.A('reupload')], None]
    else:
        return [True, ['Drag or ', html.A('select YYF file here')], None]


# In[6]:


@app.callback(
    [
    Output('grouping_upload_text', 'children'),
    Output('grouping_df', 'data'),
    ],
    [
    Input('groups_upload', 'contents'),
    ]
)
def uploadGroupingFile(grouping_file):
    if grouping_file != None:
        try:
            content_type, content_string = grouping_file.split(',')
            decoded = base64.b64decode(content_string).decode('utf-8-sig', 'ignore')
            data = [line.strip().split(',') for line in decoded.split('\n')]
            grouping_df = pd.DataFrame(data[1:], columns = data[0])
            if list(grouping_df.columns) != ['Sample','Group']:
                return [['Corrupted grouping file, ', html.A('reupload')], None]
            return [['Groups file uploaded'], grouping_df.to_dict('records')]
        except:
            return [['Corrupted groups file, ', html.A('reupload')], None]
    else:
        return [['optional: add ', html.A('groups file')], None]


# In[7]:


@app.callback(
    [
    Output('upload_status', 'children'),
    ],
    [
    Input('Analysis_start_btn', 'children'),
    Input('Analysis_start_btn', 'disabled'),
    Input('yyf_df', 'data'),
    Input('grouping_df', 'data'),
    ],
)
def updateUploadingStatus(x, y, yyf_df, grouping_df):
    if (yyf_df != None) and (grouping_df == None):
        return ['Upload a grouping file or run analysis']
    elif (yyf_df == None) and (grouping_df != None):
        return ['Upload a YYF file to run analysis']
    elif (yyf_df != None) and (grouping_df != None):
        return ['Press the button to run analysis']
    else:
        return ['Awating input files']


# In[8]:


@app.callback(
    [Output('upload_div', 'style'),
     Output('results_div', 'style'),
     Output('table', 'data'),
     Output('table', 'columns'),
     Output('table', 'row_selectable'),
     Output('referenceGenome', 'data')],
    [Input('Analysis_start_btn', 'n_clicks'),
     Input('yyf_df' , 'data')]
)
def startAnalysis(n_clicks, data):
    default_style = style = {'width' : '100%', 'display' : 'inline-block','textAlign' : 'center'}
    if data == None:
        return [default_style, {'display' : 'none'}, None, None, None, None]
    yyf_df = pd.DataFrame(data[1:], columns = data[0])
    yyf_df['ClinVar ID'] = yyf_df['ClinVar ID'].apply(lambda x : str(x).replace('.0',''))
    if n_clicks > 0 and yyf_df.empty == False:
        for col in yyf_df.columns:
            yyf_df[col] = yyf_df[col].str.replace('~~',',')
        referenceGenome = [col for col in yyf_df.columns if str(col).startswith('Coordinates - hg')][0]
        yyf_df = yyf_df[yyf_df[referenceGenome].str.len() > 3]
        shown_cols = ['Alleles', 'Coordinates - hg19', 'Coordinates - hg38', 'Gene', 'Phenotypes', 'Review Status']
        cols = [{'id': c, 'name': c} for c in yyf_df.columns if c in shown_cols]
        tooltip_data = [{column: {'value': str(value), 'type': 'markdown'} for column, value in row.items()} for row in yyf_df.to_dict('records')],
        return [{'display' : 'none'}, default_style, yyf_df.to_dict('records'), cols, 'single', referenceGenome]
    else:
        return [default_style, {'display' : 'none'}, None, None, None, None]


# In[9]:


@app.callback(
    [Output('information_window_div', 'children'),
     ],
    [Input('table', 'derived_virtual_selected_rows'),
     Input('table', 'derived_virtual_data'),
     Input('referenceGenome', 'data'),
     Input('grouping_df', 'data')]
)
def getVariantData(selected_row_index, data, referenceGenome, grouping_df):
    if selected_row_index in [[], None, ()]:
        return [[html.P('Select row for more information'), html.Img(src = 'assets/click_demo.gif')]]
    else:
        row = data[selected_row_index[0]]
        grouping_df = pd.DataFrame(grouping_df)
        alleles = row['Alleles']
        phenotypes = row['Phenotypes'].replace(',', '')
        gene = row['Gene']
        variation = row[referenceGenome]
        clinVarId = row['ClinVar ID']
        variantCols = ['Alleles', 'ClinVar ID', 'Coordinates - hg38', 'Coordinates - hg19', 'Gene', 'Phenotypes', 'Review Status']
        samples = [i for i in list(data[0].keys()) if i not in variantCols]
        informationWindow = []
        hets = [sample for sample in samples if row[sample] == '1']
        nhets = len(hets)
        homs = [sample for sample in samples if row[sample] == '2']
        nhoms = len(homs)
        
        linkStyle = {'marginRight' : 20, 'display' : 'inline-block', 'color' : "#077be2",'font-family' : 'gisha', 'fontSize' : 15, 'textAlign' : 'center' ,'borderWidth' : '2px','borderColor' : "#000044",'borderStyle' : 'groove','padding': 4,'borderRadius' : '5px'}

        if str(clinVarId).replace(' ', '').strip() != '':
            ClinVarUrl = "https://www.ncbi.nlm.nih.gov/clinvar/variation/" + str(clinVarId)
            ClinVarLink = html.P(html.A('ClinVar Page',target='_blank', href = ClinVarUrl, style = {'fontSize' : 18}), style = linkStyle)
            informationWindow.append(ClinVarLink)
        coordinates = variation.split(' ')[0]
        if 'INDEL' in variation:
            if referenceGenome.endswith('19'):
                gnomAD_Url = 'https://gnomad.broadinstitute.org/region/' + coordinates + '?dataset=gnomad_r2_1'
            elif referenceGenome.endswith('38'):
                gnomAD_Url = 'https://gnomad.broadinstitute.org/region/' + coordinates + '?dataset=gnomad_r3'
        else:
            chromosome = variation.split(':')[0]
            position = variation.split(':')[1].split(' ')[0]
            ref = variation.split(' ')[1].split('>')[0]
            alt = variation.split(' ')[1].split('>')[1]
            mutation = '-'.join([chromosome, position, ref, alt])
            if referenceGenome.endswith('19'):
                gnomAD_Url = 'https://gnomad.broadinstitute.org/variant/' + mutation + '?dataset=gnomad_r2_1'
            elif referenceGenome.endswith('38'):
                gnomAD_Url = 'https://gnomad.broadinstitute.org/variant/' + mutation + '?dataset=gnomad_r3'
        gnomADLink = html.P(html.A('gnomAD frequency',target='_blank', href = gnomAD_Url, style = {'fontSize' : 18}), style = linkStyle)
        informationWindow.append(gnomADLink)
        if grouping_df.empty == False:
            groups = [group for group in list(grouping_df['Group'].unique()) if group != None]
            groups_dict = pd.Series(grouping_df['Group'].values, index = grouping_df['Sample']).to_dict()
            counts = dict()
            for het in hets:
                try:
                    if groups_dict[het] in counts:
                        counts[groups_dict[het]] += 1
                    else:
                        counts[groups_dict[het]] = 1
                except:
                    continue
            for hom in homs:
                try:
                    if groups_dict[hom] in counts:
                        counts[groups_dict[hom]] += 1
                    else:
                        counts[groups_dict[hom]] = 1
                except:
                    continue

            bardf = pd.DataFrame(counts.items(), columns = ['Group', 'count'])
            bardf['group_total_count'] = bardf['Group'].apply(lambda x : grouping_df['Group'].tolist().count(x))
            bardf['hover_text'] = bardf['count'].astype(str) + "/" + bardf['group_total_count'].astype(str)
            bardf['percentage'] = bardf['count'] / bardf['group_total_count']
            bardf['percentage'] = bardf['percentage'].round(2) * 100
            graph = dcc.Graph(
                figure = {'data' : [
                    go.Bar(
                        x = bardf['Group'],
                        y = bardf['percentage'],
                        textposition='auto',
                        text = [str(int(round(p,2))) + "%" for p in bardf['percentage'].tolist()],
                        hovertext = ["%<br><br>" + p for p in bardf['hover_text'].tolist()],
                        opacity = 0.6
                    )
                ]}
            )
            informationWindow.append(graph)
        if nhets == 0:
            informationWindow.append(html.P('0 heterozygotes.', style = {'fontSize' : 16}))
        else:
            informationWindow.append(html.P(str(nhets) + ' heterozygotes: ' + ', '.join(hets) + '.', style = {'fontSize' : 16}))
        if nhoms == 0:
            informationWindow.append(html.P('0 homozygotes.', style = {'fontSize' : 16}))
        else:
            informationWindow.append(html.P(str(nhoms) + ' homozygotes: ' + ', '.join(homs) + '.', style = {'fontSize' : 16}))
        return [informationWindow]


# In[10]:


@app.callback(
    [Output('population_analysis_div', 'children'),
     Output('population_start_btn', 'style'),
     Output('heatmap_df', 'data')
     ],
    [Input('yyf_df', 'data'),
     Input('population_start_btn', 'n_clicks')
    ]
)
def populationAnalysis(data, n_clicks):
    if n_clicks == 0:
        return [None, {'text-transform' : 'none', 'font-size' : 24, 'margin' : 'auto'}, None]
    else:
        df = pd.DataFrame(data[1:], columns = data[0])
        df = df[df.columns[7:]]
        df.replace('','0', inplace = True)
        df = df.astype(int)
        df_euclid = pd.DataFrame(
            1 / (1 + distance_matrix(df.T, df.T)),
            columns = df.columns, index = df.columns
        )
        fig_heatmap = px.imshow(df_euclid, labels={'x': 'Sample 1', 'y':'Sample 2', 'color' : 'Similarity score'})
        fig_heatmap.update_xaxes(visible=False)
        fig_heatmap.update_yaxes(visible=False)
        heatmap = dcc.Graph(figure = fig_heatmap)
        
        explanation1 = html.P('Mouseover the heatmap to find similarity between each two samples',
                             style = {'width' : '100%', 'display' : 'inline-block','textAlign' : 'center'}
        )
        #Check most similar
        samples = [{'label': v, 'value': v} for v in df.columns]
        mostSimilarDropDown = dcc.Dropdown(options = samples,
            id = 'mostSimilarDropDown',
            placeholder = 'Show samples most similar to...',
            style = {'margin' : 'auto', 'width' : '75%', 'display' : 'inline-block','textAlign' : 'center'}
        )
        mostSimilarDiv = html.Div([html.Div(id = 'mostSimilar'), mostSimilarDropDown], style = {'width' : '100%', 'display' : 'inline-block','textAlign' : 'center'})        
        if len(samples) >= 10:
            umap_btn = html.Button('Perform UMAP clustering analysis',
                id = 'umap_start_btn', 
                n_clicks = 0,
                style = {'text-transform' : 'none', 'font-size' : 20, 'margin' : 'auto'}
            )
        else:
            umap_btn = html.Button('Perform UMAP clustering analysis',
                id = 'umap_start_btn', 
                n_clicks = 0,
                style = {'display' : 'none'}
            )
        umapDiv = html.Div([html.Div(id = 'umapResult'), umap_btn], style = {'width' : '100%', 'display' : 'inline-block','textAlign' : 'center'})        
        children = [explanation1, heatmap, html.Br(), mostSimilarDiv, html.Br(), html.Hr(), umapDiv]
        return [children, {'display' : 'none'}, df_euclid.to_dict('records')]

@app.callback(
    [Output('mostSimilar', 'children'),
     ],
    [Input('heatmap_df', 'data'),
     Input('mostSimilarDropDown', 'value')
    ]
)
def getMostSimilarSamples(data, sample):
    if sample == None:
        return [dash.no_update]
    df = pd.DataFrame(data, columns = data[0])
    df.sort_values(sample, ascending = False, inplace = True)
    df.index = df.columns
    df.drop(sample, inplace = True)
    order = list(df.index)
    title = html.H3('Samples most similar to ' + sample)
    fig = go.Figure(go.Bar(x = order, y = df[sample]))
    fig.update_layout(paper_bgcolor = 'rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')
    graph = dcc.Graph(figure = fig)
    return [html.Div([title, graph])]

@app.callback(
    [Output('umapResult', 'children'),
     Output('umap_start_btn', 'style'),
     ],
    [Input('umap_start_btn', 'n_clicks'),
     Input('yyf_df', 'data'),
    ]
)
def generateUmap(n_clicks, data):
    if n_clicks == 0:
        return [dash.no_update, dash.no_update]
    else:
        reducer = umap.UMAP(low_memory = True)
        df = pd.DataFrame(data[1:], columns = data[0])
        df = df[df.columns[7:]]
        df.replace('','0', inplace = True)
        df = df.astype(int)
        df = df.T
        data = StandardScaler().fit_transform(df)
        embedding = reducer.fit_transform(data)
        df['Sample'] = df.index
        df['UMAP 1'] = embedding[:,0]
        df['UMAP 2'] = embedding[:,1]
        df = df[['UMAP 1', 'UMAP 2', 'Sample']]
        fig_umap = px.scatter(df, x = 'UMAP 1', y = 'UMAP 2', hover_data = ['Sample'], title = 'UMAP', opacity = 0.5)
        fig_umap.update_layout(title_x = 0.5)
        Umap = dcc.Graph(figure = fig_umap)
        return [Umap, {'display' : 'none'}]


# In[ ]:


'''
import socket
s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
s.connect(('8.8.8.8', 1))
local_ip_address = s.getsockname()[0]
if __name__ == '__main__':
    app.run_server(port = 7454, host = local_ip_address)
'''
if __name__ == '__main__':
    app.run_server(debug=True)

