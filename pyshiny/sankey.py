import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

from shinywidgets import output_widget, render_widget
from typing import Callable
from shiny import Inputs, Outputs, Session, module, render, ui, reactive, req
from sklearn.decomposition import PCA
from matplotlib import cm

####################################################
#AMPCOMBI - HMM models
####################################################
@module.ui
def ampcombi_sankey_ui():
    return  ui.nav_panel(
        "Taxonomy",
       ui.p("Contig taxonomic lineage for AMP hits based on MMseqs2 classification:"),
       ui.input_selectize(
                "selected_taxa_levels",
                "Choose more than one taxonomic level:",
                choices= ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'specie'],
                selected= ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'specie'],
                multiple = True,
                ),
       #     output_widget("amp_clusters_taxa")),
        ui.card(
            output_widget("taxa_sankey")),
    )

@module.server
def ampcombi_sankey_server(
    input: Inputs,
    output: Outputs,
    session: Session,
    df: Callable[[], pd.DataFrame]
    ):
    
    @render_widget
    def taxa_sankey():
        """
        Taxonomy across the contigs with AMPs. 
        Obtained from FUNCSCAN.
        """
        if isinstance(df(), pd.DataFrame):
            return taxonomy_sankey_plot(df(),input.selected_taxa_levels())

def taxonomy_sankey_plot(data, selected_levels):
    """
    Taxonomy across the contigs with AMPs.
    This creates the SANKEY plot.
    """
    df_amp = data
    
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
            }, 
        'link': {
            }
        }
        ]
    
    # Check input presence
    if 'mmseqs_lineage_contig' in df_amp.columns:
        df_amp["lineage"] = df_amp["mmseqs_lineage_contig"].str.extract(r"(d_.*)")
        df_amp[['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'specie']] = df_amp['mmseqs_lineage_contig'].str.split(';', expand=True)
        # Remove the prefix from each column
        df_amp['kingdom'] = df_amp['kingdom'].str.replace('d_', '')
        df_amp['phylum'] = df_amp['phylum'].str.replace('p_', '')
        df_amp['class'] = df_amp['class'].str.replace('c_', '')
        df_amp['order'] = df_amp['order'].str.replace('o_', '')
        df_amp['family'] = df_amp['family'].str.replace('f_', '')
        df_amp['genus'] = df_amp['genus'].str.replace('g_', '')
        df_amp['specie'] = df_amp['specie'].str.replace('s_', '')
        # Remove the letters used in GTDB formating
        df_amp['phylum'] = df_amp['phylum'].str.replace(r'\s[A-Z](?!\w)', '', regex=True)
        df_amp['class'] = df_amp['class'].str.replace(r'\s[A-Z](?!\w)', '', regex=True)
        df_amp['order'] = df_amp['order'].str.replace(r'\s[A-Z](?!\w)', '', regex=True)
        df_amp['family'] = df_amp['family'].str.replace(r'\s[A-Z](?!\w)', '', regex=True)
        df_amp['genus'] = df_amp['genus'].str.replace(r'\s[A-Z](?!\w)', '', regex=True)
        df_amp['specie'] = df_amp['specie'].str.replace(r'\s[A-Z](?!\w)', '', regex=True)
        # Remove the genus from specie column
        df_amp['specie_only'] = df_amp['specie'].str.split(' ', n=1).str[1]
    else:
        return None
    
    # Start plotting:
    # Remove NANs and roots in kingdom
    unknown_contigs = df_amp[~df_amp['kingdom'].isin(['Bacteria', 'Archaea'])]
    df_amp = df_amp[df_amp['kingdom'].isin(['Bacteria', 'Archaea'])]
    # Grab the first 100 genus and filter the df to it :: this solves the complexity of the sankey plot - change here the number of genera to illustrate
    genus_counts = df_amp['genus'].value_counts().sort_values(ascending=False)
    top_abundant = genus_counts.index[:2].tolist()
    df_amp = df_amp[df_amp['genus'].isin(top_abundant)]
    # Remove alphabets from GTDB db in species name to reduce heterogneity
    df_amp['specie_only'] = df_amp['specie_only'].str.split(' ', n=1).str[0]
    df_amp

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
    df_amp_GS['specie'] = df_amp['specie_only']
    # Count the number of times the combinations are there for specie+genus
    df_amp_GS = df_amp_GS.groupby(['genus', 'specie'], sort=False, dropna=False).size().reset_index(name='count')
    # Replace string values with unique numbers in new columns
    df_amp_GS['genus_id'] = pd.factorize(df_amp_GS['genus'])[0] + max_value
    max_value = (df_amp_GS['genus_id'].max()+1) # grab the maximum value from the column
    df_amp_GS['specie_id'] = pd.factorize(df_amp_GS['specie'])[0] + max_value
    # Stack them vertically
    df_id = pd.melt(df_amp_GS, value_vars=['genus_id', 'specie_id'], var_name='Category', value_name='node_id')
    df_name = pd.melt(df_amp_GS, value_vars=['genus', 'specie'], var_name='Category', value_name='label')
    # Merge the two melted DataFrames based on the index
    df_amp_nodes_GS = pd.concat([df_name, df_id['node_id']], axis=1)
    # Drop duplicates
    df_amp_nodes_GS.drop_duplicates(keep='first', inplace=True, ignore_index=True)

    ########################
    # (7) NODES labels
    ########################
    df_amp_nodes = pd.DataFrame()
    # Create a list of DataFrames to concatenate all nodes labels
    df_list = []

    # Conditional inclusion based on chosen taxonomy levels
    if 'kingdom' in selected_levels and 'phylum' in selected_levels:
        df_list.append(df_amp_nodes_KP)
    if 'phylum' in selected_levels and 'class' in selected_levels:
        df_list.append(df_amp_nodes_PC)
    if 'class' in selected_levels and 'order' in selected_levels:
        df_list.append(df_amp_nodes_CO)
    if 'order' in selected_levels and 'family' in selected_levels:
        df_list.append(df_amp_nodes_OF)
    if 'family' in selected_levels and 'genus' in selected_levels:
        df_list.append(df_amp_nodes_FG)
    if 'genus' in selected_levels and 'specie' in selected_levels:
        df_list.append(df_amp_nodes_GS)

    # Concatenate only the selected DataFrames
    df_amp_nodes = pd.concat(df_list, ignore_index=True)
    
    # Drop unnecessary columns
    df_amp_nodes.drop_duplicates(keep='first', inplace=True, ignore_index=True)
    values_to_remove = [ 'Category', 'node_id'] 
    df_amp_nodes.drop(values_to_remove,axis=1, inplace=True)
    # Remove empty rows
    df_amp_nodes = df_amp_nodes.dropna()
    # Assign a number from 0 to inf
    df_amp_nodes['node_id'] = range(len(df_amp_nodes))
    # Convert column to list 
    column_list = df_amp_nodes['label'].tolist()
    # Replace the 'label' item in the 'node' dictionary with the column list
    data[0]['node']['label'] = column_list
    # Generate distinct colors according to node_no
    nodes_no = len(column_list) #number of nodes
    # Colour palette
    palette = cm.get_cmap('viridis_r', nodes_no)
    color_list = [
        f"rgba({int(r * 255)}, {int(g * 255)}, {int(b * 255)}, 0.8)"
        for r, g, b, _ in palette(np.linspace(0, 1, nodes_no))
    ]
    data[0]['node']['color'] = color_list

    ########################
    # (8) COUNTS labels
    ########################
    # replace the contents of xx_id with the correct digits
    # Remove '_id' columns if present in each valid DataFrame : drop unnecessary columns
    dataframes = [df_amp_KP, df_amp_PC, df_amp_CO, df_amp_OF, df_amp_FG, df_amp_GS]
    for df in dataframes:
        if df is not None: 
            columns_to_remove = [col for col in df.columns if '_id' in col]
            df.drop(columns_to_remove, axis=1, inplace=True)
    # Rename columns only for valid DataFrames : stack and rename the columns
    for df in dataframes:
        if df is not None:
            df.columns = ['source', 'target', 'count']
        
    # Conditional inclusion based on chosen taxonomy levels
    if 'kingdom' in selected_levels and 'phylum' in selected_levels:
        df_list.append(df_amp_KP)
    if 'phylum' in selected_levels and 'class' in selected_levels:
        df_list.append(df_amp_PC)
    if 'class' in selected_levels and 'order' in selected_levels:
        df_list.append(df_amp_CO)
    if 'order' in selected_levels and 'family' in selected_levels:
        df_list.append(df_amp_OF)
    if 'family' in selected_levels and 'genus' in selected_levels:
        df_list.append(df_amp_FG)
    if 'genus' in selected_levels and 'specie' in selected_levels:
        df_list.append(df_amp_GS)

    # Concatenate only the selected DataFrames
    df_amp_counts = pd.concat(df_list,  axis=0,ignore_index=True)
    
    # Remove rows that have no values in two columns
    df_amp_counts = df_amp_counts.dropna(subset=['source', 'target'])
    # Labels and digits to dict
    label_dict = df_amp_nodes.set_index('label')['node_id'].to_dict()
    # Replace the values in the 'source' column with digits
    df_amp_counts['source_digit'] = df_amp_counts['source'].replace(label_dict).astype(int)
    # Replace the values in the 'target' column with digits
    df_amp_counts['target_digit'] = df_amp_counts['target'].map(label_dict).astype(int)
    # Convert columns to list 
    column_list = df_amp_counts['count'].tolist()
    data[0]['link']['value'] = column_list
    column_list = df_amp_counts['target_digit'].tolist()
    data[0]['link']['target'] = column_list
    column_list = df_amp_counts['source_digit'].tolist()
    data[0]['link']['source'] = column_list
    # Make sure the values in the list are integers 
    data[0]['link']['source'] = [int(src) for src in data[0]['link']['source']]
    data[0]['link']['target'] = [int(src) for src in data[0]['link']['target']]
    data[0]['link']['value'] = [int(src) for src in data[0]['link']['value']]
    
    ########################
    # (9) Prepare for the Sankey diagram
    ########################
    # Adapted from https://python-graph-gallery.com/sankey-diagram-with-python-and-plotly/
    opacity = 0.4

    # Change 'magenta' to its 'rgba' value to add opacity
    data[0]['link']['color'] = [data[0]['node']['color'][src].replace("0.8", str(opacity))
                                        for src in data[0]['link']['target']]
    # Change opacity of the links
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

    return fig.update_layout(
    font=dict(family="Arial", size=14, color="black"),
    hoverlabel=dict(font=dict(family='Arial', size=14))
    )
