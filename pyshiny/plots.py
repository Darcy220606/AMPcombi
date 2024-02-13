from pandas import DataFrame
import pandas as pd
from pathlib import Path
from shiny.types import ImgData
import matplotlib.pyplot as plt
import matplotlib.style as style
from upsetplot import UpSet

import plotly.express as px
import plotly.graph_objects as go

import urllib, json
import urllib.request
import json
import requests
import pandas as pd
import colorsys
from IPython.display import display, HTML

import py3Dmol

def plot_ampcombi_upset(df: DataFrame):
    """
    This is activated by the user. The png/html/pdf corresponding to the
    tool comparison plot (UPset plot) for ampcombi.        
    """
    selected_columns = ['prob_ampir',
                        'prob_neubi',
                        'HMM_model',
                        'prob_ampgram',
                        'prob_amptransformer',
                        'prob_macrel',
                        'prob_amplify',
                        'CDS_id']

    df = df[selected_columns].copy()

    # convert integers to boolean
    columns_to_check = ['prob_ampir', 
     'prob_neubi',
     'prob_ampgram',
     'prob_amptransformer',
     'prob_macrel',
     'prob_amplify']
    df[columns_to_check] = df[columns_to_check].astype(float)
    #for column in columns_to_check:
        #df[column] = df[column].map(lambda x: True if x > 0.0 and pd.notna(x) else False)
        #df[column] = df[column].map(lambda x: False if pd.isna(x) or x == '' or x == 0 or x == 0.0 else True)
        #df[column] = df[column].map(lambda x: True if x > 0.0 else False)
    for column in columns_to_check:
        df[column] = df[column].map(lambda x: True if pd.notna(x) and x != '' and float(x) > 0.0 else False)
    
    # for HMM_model seperate due to category and 0.0
    #df['HMM_model'] = df['HMM_model'].map(lambda x: False if pd.isnull(x) or x == '0' or x == '0.0' else True)
    df['HMM_model'] = df['HMM_model'].map(lambda x: False if pd.isnull(x) or x in ['0', '', '0.0'] else True)
    df = pd.DataFrame(df)

    #rename the boolean cols
    rename_dict = {'prob_ampir' : 'ampir', 
                   'prob_neubi' : 'neubi',
                   'prob_ampgram' : 'ampgram',
                   'prob_amptransformer' : 'amptransformer',
                   'prob_macrel' : 'macrel',
                   'prob_amplify' : 'amplify',
                   'HMM_model':'HMMsearch'}
    df = df.rename(columns=rename_dict) 
    # set the columns as categories
    df.columns = pd.CategoricalIndex(df.columns, categories=['ampir', 'neubi', 'HMMsearch', 'amptransformer', 'ampgram', 'macrel', 'amplify', 'CDS_id'], ordered=True)

    # convert the DataFrame to a Series with multi-level index
    series = df.set_index(['ampir', 'neubi', 'HMMsearch', 'amptransformer', 'ampgram', 'macrel', 'amplify'])['CDS_id'].astype(str).sort_index().squeeze()

    # create the UpSet plot
    upset_data = series.reset_index().rename(columns={0: 'CDS_id'})
    levels = ['ampir', 'neubi', 'HMMsearch', 'amptransformer', 'ampgram', 'macrel', 'amplify']
    upset_data = upset_data.set_index(levels)['CDS_id']
    
    # plot
    # for dark shade : ohne dark background justremove the 'with'
    with plt.style.context('dark_background'):
    #style.use('default') 
        upset = UpSet(upset_data, subset_size='count', show_counts=True, 
                  facecolor="lightblue", #lightblue
                  other_dots_color=.4,
             shading_color=.5, totals_plot_elements=5,element_size=30) #shading_color="lightgray",element_size=30

        fig = plt.figure(figsize=(8, 6))
        plot = upset.plot()

        # Save the plot 
        #plt.savefig('ampcombi_upset_plot.pdf')
        #plt.savefig('ampcombi_upset_plot.png',dpi=300)
    return plt.show()

def plot_ampcombi_cluster(df: DataFrame):
    # Count the number of occurrences of each unique cluster_id
    unique_counts = df['cluster_id'].value_counts().reset_index()

    # Rename the columns to 'cluster_id' and 'counts'
    unique_counts = unique_counts.rename(columns={'index': 'cluster_id', 'cluster_id': 'counts'})
    unique_counts['cluster_id'] = unique_counts['cluster_id'].astype(str)

    scatterplot = px.bar(
        data_frame=unique_counts,x="cluster_id", y="counts", text="counts").update_layout(
        title={"text": "Clustered AMP Hits"},
        yaxis_title="Number of hits in cluster",
        xaxis_title="Cluster ID",
    )
    scatterplot.update_traces(textposition='auto') #write on bars
    scatterplot.update_traces(marker_color='darkgray') #colours the bars
    return scatterplot

def plot_ampcombi_sankey(filtered):
    # set up data and layout required by sankey plots
    data = [
        {
        'type': 'sankey', 
        'domain': {'x': [0, 1], 'y': [0, 1]}, 
        'orientation': 'h', 
        'valueformat': '.0f', 
        'valuesuffix': 'TWh', 
        'node': {
            'pad': 15, 
            'thickness': 15, 
            'line': {'color': 'black', 'width': 0.5}, 
            #'label': ["Agricultural 'waste'", 'Bio-conversion', 'Liquid', 'Losses', 'Solid', 'Gas', 'Biofuel imports', 'Biomass imports', 'Coal imports', 'Coal', 'Coal reserves', 'District heating', 'Industry', 'Heating and cooling - commercial', 'Heating and cooling - homes', 'Electricity grid', 'Over generation / exports', 'H2 conversion', 'Road transport', 'Agriculture', 'Rail transport', 'Lighting & appliances - commercial', 'Lighting & appliances - homes', 'Gas imports', 'Ngas', 'Gas reserves', 'Thermal generation', 'Geothermal', 'H2', 'Hydro', 'International shipping', 'Domestic aviation', 'International aviation', 'National navigation', 'Marine algae', 'Nuclear', 'Oil imports', 'Oil', 'Oil reserves', 'Other waste', 'Pumped heat', 'Solar PV', 'Solar Thermal', 'Solar', 'Tidal', 'UK land based bioenergy', 'Wave', 'Wind'], 
            #'color': ['rgba(31, 119, 180, 0.8)', 'rgba(255, 127, 14, 0.8)', 'rgba(44, 160, 44, 0.8)', 'rgba(214, 39, 40, 0.8)', 'rgba(148, 103, 189, 0.8)', 'rgba(140, 86, 75, 0.8)', 'rgba(227, 119, 194, 0.8)', 'rgba(127, 127, 127, 0.8)', 'rgba(188, 189, 34, 0.8)', 'rgba(23, 190, 207, 0.8)', 'rgba(31, 119, 180, 0.8)', 'rgba(255, 127, 14, 0.8)', 'rgba(44, 160, 44, 0.8)', 'rgba(214, 39, 40, 0.8)', 'rgba(148, 103, 189, 0.8)', 'rgba(140, 86, 75, 0.8)', 'rgba(227, 119, 194, 0.8)', 'rgba(127, 127, 127, 0.8)', 'rgba(188, 189, 34, 0.8)', 'rgba(23, 190, 207, 0.8)', 'rgba(31, 119, 180, 0.8)', 'rgba(255, 127, 14, 0.8)', 'rgba(44, 160, 44, 0.8)', 'rgba(214, 39, 40, 0.8)', 'rgba(148, 103, 189, 0.8)', 'rgba(140, 86, 75, 0.8)', 'rgba(227, 119, 194, 0.8)', 'rgba(127, 127, 127, 0.8)', 'rgba(188, 189, 34, 0.8)', 'rgba(23, 190, 207, 0.8)', 'rgba(31, 119, 180, 0.8)', 'rgba(255, 127, 14, 0.8)', 'rgba(44, 160, 44, 0.8)', 'rgba(214, 39, 40, 0.8)', 'rgba(148, 103, 189, 0.8)', 'magenta', 'rgba(227, 119, 194, 0.8)', 'rgba(127, 127, 127, 0.8)', 'rgba(188, 189, 34, 0.8)', 'rgba(23, 190, 207, 0.8)', 'rgba(31, 119, 180, 0.8)', 'rgba(255, 127, 14, 0.8)', 'rgba(44, 160, 44, 0.8)', 'rgba(214, 39, 40, 0.8)', 'rgba(148, 103, 189, 0.8)', 'rgba(140, 86, 75, 0.8)', 'rgba(227, 119, 194, 0.8)', 'rgba(127, 127, 127, 0.8)']
            }, 
        'link': {
            #'source': [0, 1, 1, 1, 1, 6, 7, 8, 10, 9, 11, 11, 11, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 23, 25, 5, 5, 5, 5, 5, 27, 17, 17, 28, 29, 2, 2, 2, 2, 2, 2, 2, 2, 34, 24, 35, 35, 36, 38, 37, 39, 39, 40, 40, 41, 42, 43, 43, 4, 4, 4, 26, 26, 26, 44, 45, 46, 47, 35, 35], 
            #'target': [1, 2, 3, 4, 5, 2, 4, 9, 9, 4, 12, 13, 14, 16, 14, 17, 12, 18, 19, 13, 3, 20, 21, 22, 24, 24, 13, 3, 26, 19, 12, 15, 28, 3, 18, 15, 12, 30, 18, 31, 32, 19, 33, 20, 1, 5, 26, 26, 37, 37, 2, 4, 1, 14, 13, 15, 14, 42, 41, 19, 26, 12, 15, 3, 11, 15, 1, 15, 15, 26, 26], 
            #'value': [124.729, 0.597, 26.862, 280.322, 81.144, 35, 35, 11.606, 63.965, 75.571, 10.639, 22.505, 46.184, 104.453, 113.726, 27.14, 342.165, 37.797, 4.412, 40.858, 56.691, 7.863, 90.008, 93.494, 40.719, 82.233, 0.129, 1.401, 151.891, 2.096, 48.58, 7.013, 20.897, 6.242, 20.897, 6.995, 121.066, 128.69, 135.835, 14.458, 206.267, 3.64, 33.218, 4.413, 14.375, 122.952, 500, 139.978, 504.287, 107.703, 611.99, 56.587, 77.81, 193.026, 70.672, 59.901, 19.263, 19.263, 59.901, 0.882, 400.12, 46.477, 525.531, 787.129, 79.329, 9.452, 182.01, 19.013, 289.366, 100, 100] 
            #'color': ['rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(33,102,172,0.35)', 'rgba(178,24,43,0.35)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'rgba(0,0,96,0.2)', 'lightgreen', 'goldenrod'] 
            #'label': ['stream 1', '', '', '', 'stream 1', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'stream 1', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'Old generation plant (made-up)', 'New generation plant (made-up)', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
            }
        }
        ]

    df_amp = filtered
    # grab the cluster id if set by the user 
    #if input.clusters_id_tax() is not None:
    #    df_amp = df_amp[df_amp['cluster_id'] == input.clusters_id_tax()]

    # grab the contig GTDB classifications
    df_amp[['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'specie']] = df_amp['mmseqs_lineage_contig'].str.split(';', expand=True)
    # remove the prefix from each column
    df_amp['kingdom'] = df_amp['kingdom'].str.replace('d_', '')
    df_amp['phylum'] = df_amp['phylum'].str.replace('p_', '')
    df_amp['class'] = df_amp['class'].str.replace('c_', '')
    df_amp['order'] = df_amp['order'].str.replace('o_', '')
    df_amp['family'] = df_amp['family'].str.replace('f_', '')
    df_amp['genus'] = df_amp['genus'].str.replace('g_', '')
    df_amp['specie'] = df_amp['specie'].str.replace('s_', '')
    # remove the letters used in GTDB formating
    df_amp['phylum'] = df_amp['phylum'].str.replace(r'\s[A-Z](?!\w)', '', regex=True)
    df_amp['class'] = df_amp['class'].str.replace(r'\s[A-Z](?!\w)', '', regex=True)
    df_amp['order'] = df_amp['order'].str.replace(r'\s[A-Z](?!\w)', '', regex=True)
    df_amp['family'] = df_amp['family'].str.replace(r'\s[A-Z](?!\w)', '', regex=True)
    df_amp['genus'] = df_amp['genus'].str.replace(r'\s[A-Z](?!\w)', '', regex=True)
    df_amp['specie'] = df_amp['specie'].str.replace(r'\s[A-Z](?!\w)', '', regex=True)
    # remove the genus from specie column
    df_amp['specie_mod'] = df_amp['specie'].str.split(' ', n=1).str[1]
    
    ########################
    # (1-5) Kingdom//Phylum//Class/Genus
    ########################
    def create_taxonomy_nodes(df, level1, level2):
        """
        Create a dataframe of taxonomy nodes for the given levels:
        Kingdom/Phylum/Class/Order/Family/Genus/Specie.
        Parameters/inputs required:
            df (pandas.DataFrame): the input dataframe containing the taxonomy data.
            level1 (str): the name of the first level column.
            level2 (str): the name of the second level column.
        """
        df_tax = pd.DataFrame()
        df_tax[level1] = df[level1]
        df_tax[level2] = df[level2]
        df_tax = df_tax.groupby([level1, level2], sort=False, dropna=False).size().reset_index(name='count')
        df_tax[f"{level1}_id"] = pd.factorize(df_tax[level1])[0]
        max_value = (df_tax[f"{level1}_id"].max()+1)
        df_tax[f"{level2}_id"] = pd.factorize(df_tax[level2])[0] + max_value
        df_id = pd.melt(df_tax, value_vars=[f"{level1}_id", f"{level2}_id"], var_name='Category', value_name='node_id')
        df_name = pd.melt(df_tax, value_vars=[level1, level2], var_name='Category', value_name='label')
        df_tax_nodes = pd.concat([df_name, df_id['node_id']], axis=1)
        df_tax_nodes.drop_duplicates(keep='first', inplace=True, ignore_index=True)
        return df_tax, df_tax_nodes, max_value
    df_amp_KP, df_amp_nodes_KP, max_value = create_taxonomy_nodes(df_amp, 'kingdom', 'phylum')
    df_amp_PC, df_amp_nodes_PC, max_value = create_taxonomy_nodes(df_amp, 'phylum', 'class')
    df_amp_CO, df_amp_nodes_CO, max_value = create_taxonomy_nodes(df_amp, 'class', 'order')
    df_amp_OF, df_amp_nodes_OF, max_value = create_taxonomy_nodes(df_amp, 'order', 'family')
    df_amp_FG, df_amp_nodes_FG, max_value = create_taxonomy_nodes(df_amp, 'family', 'genus')
    
    ########################
    # (6) Genus // Specie (sp.)
    ########################
    df_amp_GS = pd.DataFrame()
    df_amp_GS['genus'] = df_amp['genus']
    df_amp_GS['specie'] = df_amp['specie_mod']
    # count the number of times the combinations are there for specie+genus
    df_amp_GS = df_amp_GS.groupby(['genus', 'specie'], sort=False, dropna=False).size().reset_index(name='count')
    # replace string values with unique numbers in new columns
    df_amp_GS['genus_id'] = pd.factorize(df_amp_GS['genus'])[0] + max_value
    max_value = (df_amp_GS['genus_id'].max()+1) # grab the maximum value from the column
    df_amp_GS['specie_id'] = pd.factorize(df_amp_GS['specie'])[0] + max_value
    # stack them vertically
    df_id = pd.melt(df_amp_GS, value_vars=['genus_id', 'specie_id'], var_name='Category', value_name='node_id')
    df_name = pd.melt(df_amp_GS, value_vars=['genus', 'specie'], var_name='Category', value_name='label')
    # merge the two melted DataFrames based on the index
    df_amp_nodes_GS = pd.concat([df_name, df_id['node_id']], axis=1)
    # drop duplicates
    df_amp_nodes_GS.drop_duplicates(keep='first', inplace=True, ignore_index=True)
    
    ########################
    # (7) NODES labels
    ########################
    df_amp_nodes = pd.DataFrame()
    # concatenate all nodes' labels
    df_amp_nodes = pd.concat([df_amp_nodes_KP, df_amp_nodes_PC,
                              df_amp_nodes_CO,df_amp_nodes_OF,
                              df_amp_nodes_FG,df_amp_nodes_GS])
    # drop unnecessary columns
    df_amp_nodes.drop_duplicates(keep='first', inplace=True, ignore_index=True)
    values_to_remove = [ 'Category', 'node_id'] 
    df_amp_nodes.drop(values_to_remove,axis=1, inplace=True)
    # remove empty rows
    df_amp_nodes = df_amp_nodes.dropna()
    # assign a number from 0 to inf
    df_amp_nodes['node_id'] = range(len(df_amp_nodes))
    # replace data[0]['node']['label'] / ['color'] with user input
    # convert column to list 
    column_list = df_amp_nodes['label'].tolist()
    # replace the 'label' item in the 'node' dictionary with the column list
    data[0]['node']['label'] = column_list
    # generate distinct colors according to node_no
    nodes_no = len(column_list) #number of nodes
    color_list = []
    for i in range(nodes_no):
        hue = i / nodes_no
        rgb = colorsys.hsv_to_rgb(hue, 1.0, 1.0)
        rgba = tuple(int(255 * x) for x in rgb) + (0.8,)
        color_string = f"rgba{rgba}"
        color_list.append(color_string)
    # replace the 'color' item in the 'node' dictionary with the color list
    # color_list =  ['blue'] * 137
    data[0]['node']['color'] = color_list
    
    ########################
    # (8) COUNTS labels
    ########################
    # replace the contents of xx_id with the correct digits
    # drop unnecessary columns
    dataframes = [df_amp_KP, df_amp_PC, df_amp_CO, df_amp_OF, df_amp_FG, df_amp_GS]
    for df in dataframes:
        columns_to_remove = [col for col in df.columns if '_id' in col]
        df.drop(columns_to_remove, axis=1, inplace=True)
    # stack columns in order
    # rename columns
    dataframes = [df_amp_KP, df_amp_PC, df_amp_CO, df_amp_OF, df_amp_FG, df_amp_GS]
    for df in dataframes:
        df.columns = ['source', 'target', 'count']
    df_amp_counts = pd.concat([df_amp_KP, df_amp_PC, df_amp_CO,
                               df_amp_OF, df_amp_FG, df_amp_GS], axis=0, ignore_index=True)
    # remove rows that have no values in two columns
    df_amp_counts = df_amp_counts.dropna(subset=['source', 'target'])
    # labels and digits to dict
    label_dict = df_amp_nodes.set_index('label')['node_id'].to_dict()
    # replace the values in the 'source' column with digits
    df_amp_counts['source_digit'] = df_amp_counts['source'].replace(label_dict).astype(int)
    # replace the values in the 'target' column with digits
    #df_amp_counts['target_digit'] = df_amp_counts['target'].replace(label_dict)
    df_amp_counts['target_digit'] = df_amp_counts['target'].map(label_dict).astype(int)
    # replace data[0]['link']['source'] / ['target'] / ['value'] / ['color'] with user input
    # convert columns to list 
    column_list = df_amp_counts['count'].tolist()
    data[0]['link']['value'] = column_list
    column_list = df_amp_counts['target_digit'].tolist()
    data[0]['link']['target'] = column_list
    column_list = df_amp_counts['source_digit'].tolist()
    data[0]['link']['source'] = column_list
    # make sure the values in the list are integers 
    data[0]['link']['source'] = [int(src) for src in data[0]['link']['source']]
    data[0]['link']['target'] = [int(src) for src in data[0]['link']['target']]
    data[0]['link']['value'] = [int(src) for src in data[0]['link']['value']]
    
    ########################
    # (9) Prepare for the Sankey diagram
    ########################
    # Adapted from https://python-graph-gallery.com/sankey-diagram-with-python-and-plotly/
    # override gray link colors with 'source' colors
    opacity = 0.4
    # change 'magenta' to its 'rgba' value to add opacity
    data[0]['node']['color'] = ['rgba(255,0,255, 0.8)' if color == "magenta" else color for color in data[0]['node']['color']]
    data[0]['link']['color'] = [data[0]['node']['color'][src].replace("0.8", str(opacity))
                                        for src in data[0]['link']['source']]
    # change opacity of the links
    data[0]['link']['color'] = [color.replace('0.4', '0.2') for color in data[0]['link']['color']]
    
    fig = go.Figure(data=[go.Sankey(
        valueformat = ".0f",
        valuesuffix = "TWh",
        # Define nodes
        node = dict(
          pad = 15,
          thickness = 15,
          line = dict(color = "black", width = 0.5),
          label =  data[0]['node']['label'],
          color =  data[0]['node']['color'],
          hovertemplate='%{label}<extra></extra>',
          hoverlabel=dict(font=dict(family='Arial')),
          ),
        # Add links
        link = dict(
          source =  data[0]['link']['source'],
          target =  data[0]['link']['target'],
          value =  data[0]['link']['value'],
          label =  data[0]['link']['value'],
          color =  data[0]['link']['color']
    ))],)             
    fig.update_layout(title_text="Contig taxonomic lineage for AMP hits based on MMseqs2 classification - Sankey plot ",
                      font_size=12)
    #show plot:
    #fig.show()

    fig.update_layout(
        font=dict(family="Arial", size=5, color="black"),
        hoverlabel=dict(font=dict(family='Arial'))
    )
    #save the plot in html format :
    #fig.write_html("ampcombi_sankey_plot.html")
    #fig.write_image("ampcombi_sankey_plot.pdf")
    #fig.write_image("ampcombi_sankey_plot.png", scale=5)
    return fig.show() 

def plot_ampcombi_pdb(pdb_files):
    # define a color map
    color_map = ['red', 'blue', 'green', 'yellow', 'orange', 'purple']

    # initialize the viewer
    view = py3Dmol.view(width=600, height=500)

    # create a legend
    legend_html = '<div style="position: absolute; top: 10px; left: 10px; background-color: white; padding: 10px; border: 1px solid black; border-radius: 5px;">'
    for idx, pdb_file in enumerate(pdb_files):
        legend_html += f'<div><span style="color: {color_map[idx % len(color_map)]};">■</span> {pdb_file}</div>'
    legend_html += '</div>'
    # create overlaid structures
    for idx, pdb_file in enumerate(pdb_files):
        with open(pdb_file) as ifile:
            pdb_data = ifile.read()
        view.addModel(pdb_data, 'pdb')
        view.setStyle({'model': idx}, {'cartoon': {'color': color_map[idx % len(color_map)]}})

    # zoom to fit all models
    view.zoomTo()

    # set background color
    view.setBackgroundColor('black')

    # write the structure
    output_file = 'ampcombi_pdb_structures.html'
    view.write_html(output_file)
    # add the legend HTML to the generated HTML file
    with open(output_file, 'a') as ofile:
        ofile.write(legend_html) 
    
    # render the HTML file in the default browser
    import webbrowser
    webbrowser.open_new_tab(output_file)
    #webbrowser.open(output_file)
