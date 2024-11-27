import pandas as pd
import numpy as np

from typing import Callable
from shiny import Inputs, Outputs, Session, module, render, ui, reactive, req
from shinywidgets import output_widget, render_widget, register_widget

#import io

####################################################
#AMPCOMBI - main table
####################################################
@module.ui
def ampcombi_table_ui():
    return  ui.nav_panel(
        "Table", # Name of tab
        # download rows selected : table tab
        ui.row(
            ui.card(
                ui.p("Download only the selected rows:"),
                ui.download_button("download_ampcombi_table_rows", "Download 'ampcombi_selected_rows.tsv'", class_="btn btn-info"),
                ),
            ),
        ui.row(
            ui.card(
                ui.p("Download the current filtered table:"),
                ui.download_button("download_filtered_data", "Download 'ampcombi_filtered_table.tsv'", class_="btn btn-info"),
                ),
            ),
        ui.p("Rows selected:", style="font-size: 10px;"),
        ui.output_text("ampcombi_table_rows", inline=True),
        ui.output_data_frame("ampcombi_table_dataframe"),
        #output_widget("data_grid")        
    )
    
@module.server
def ampcombi_table_server(
    input: Inputs,
    output: Outputs,
    session: Session,
    df: Callable[[], pd.DataFrame],
    ):
    selected_rows = reactive.Value([])
    
    @output
    @render.data_frame
    def ampcombi_table_dataframe():
        """"
        AMPCOMBI: render dataframe in a table
        """
        if isinstance(df(), pd.DataFrame):
            # Render grid table
            data_grid = render.DataGrid(df(),                                           
                                            row_selection_mode= 'multiple',
                                            width="100%",
                                            height="1000px",
                                            filters=True,
                                            editable=True) # This is UI-only filtering
            return data_grid
        else: None    
    
    @output
    @render.text
    def ampcombi_table_rows():
        """
        AMPCOMBI: prints the row numbers selected by user
        """
        selected = input.ampcombi_table_dataframe_selected_rows()
        # Reorder
        sorted_selected = sorted(selected)
        sorted_l = ', '.join(str(i) for i in sorted_selected)
        return sorted_l

    @output
    @render.download(filename=lambda: "ampcombi_table_selected_rows.tsv")
    async def download_ampcombi_table_rows():
        indices = list(input.ampcombi_table_dataframe_selected_rows())
        selected_rows = df().iloc[indices]
        yield selected_rows.to_csv('ampcombi_selected_rows.tsv', sep='\t', index=False)

    @reactive.calc
    def filtered_df():
        """
        Access the data view of the rendered data frame.
        """
        data_selected = ampcombi_table_dataframe.data_view(selected=False)
        req(not data_selected.empty)
        return data_selected
    
    @output
    @render.download(filename=lambda: "ampcombi_filtered_table.tsv")
    async def download_filtered_data():
        """
        Allow the user to download the filtered table as a TSV file.
        """
        yield filtered_df().to_csv('ampcombi_filtered_table.tsv', sep='\t', index=False) 
    
def datatable_filter(data, amp_tools, prediction_score, amp_lengths, min_tool_n, remove_signal_peptide):
    """
    Filters the given datatable, applying column renaming, type conversion, 
    and tool-based filtering.
    """
    # Rename columns for tools
    rename_columns = {
     "amp_length":"AMP Length",
     "molecular_weight":"MwT",
     "isoelectric_point":"IP",
     "hydrophobicity":"Hydrophobicity",
     "helix_fraction":"Helix fraction",
     "turn_fraction":"Turn fraction",
     "sheet_fraction":"Sheet fraction",
     "prob_ampir":"AMPir",
     "prob_ampgram":"AMPgram",
     "prob_neubi":"NEUBI",
     "prob_amptransformer":"AMPtransformer",
     "prob_amplify":"AMPlify",
     "prob_macrel":"MACREL",
     "evalue_hmmer":"HMMsearch",
     }
    data.rename(columns=rename_columns, inplace=True)

    # Combine all three identifiers if present
    columns_to_check = ['sample_id', 'contig_id', 'CDS_id']
    if all(col in data.columns for col in columns_to_check):
        data['seq_header_full'] = data['sample_id'] + '!' + data['contig_id'] + '!' + data['CDS_id']
    # Drop seq_headers column if it exists
    if 'seq_headers' in data.columns:
        data.drop(columns=['seq_headers'], inplace=True)
    # Change data type for 'evalue_hmmer' and format it
    if 'HMMsearch' in data.columns:
        data['HMMsearch'] = data['HMMsearch'].astype(float).apply(lambda x: '{:.2e}'.format(x))
    # Convert 'HMM_model' to categorical
    if 'HMM_model' in data.columns:
        data['HMM_model'] = data['HMM_model'].astype('category')
    # Rename 'cluster_id' and change its type
    if 'cluster_id' in data.columns:
        data['cluster_id'] = data['cluster_id'].astype(str)
        # data.rename(columns={'cluster_id': 'cluster_ID'}, inplace=True)
    
    # Add 'User specific details' column if not present
    if 'User specific details' not in data.columns:
        data['User specific details'] = pd.NA

    # Filter based on amp lengths
    if amp_lengths != 1500:
        data['AMP Length'] = data['aa_sequence'].apply(len)
        data = data[data['AMP Length'] < amp_lengths]
    else: 
        data['AMP Length'] = data['aa_sequence'].apply(len)

    # Filter based on prediction score
    if prediction_score != 0:
        filtered_dfs_2 = []
        # Iterate over each tool in amp_tools
        for tool in amp_tools:
            # Skip the 'HMM_model' tool if it is in the list
            if tool == 'HMM_model':
                continue
            if tool in data.columns:
                filtered_df_2 = data[data[tool] >= prediction_score]
                filtered_dfs_2.append(filtered_df_2)
        # Concatenate all filtered dataframes back together
        if filtered_dfs_2:
            data = pd.concat(filtered_dfs_2).drop_duplicates().reset_index(drop=True)

    # Filter based on min tool amp within clusters (IF clusters PRESENT)
    if 'cluster_id' in data.columns:
        data = amps_in_more_than_x_tools(data, amp_tools, min_tool_n)
           
    # Filter based on signal peptide within cluster numbers (IF clusters and signal peptide PRESENT)
    if ((all(col in data.columns for col in ['cluster_id', 'signal_peptide'])) and (remove_signal_peptide == True) ) :
        data = amps_with_signalpeptide(data)
        #data = signal_peptide_in_more_than_x_clusters(data, min_signal_peptide_n)
        
        
    # Return original dataframe if no tools are selected
    if not amp_tools:
        return data
    
    # Filter based on selected tools
    filtered_tables = []
    for tool in amp_tools:
        filtered_data = data[(data[tool] != 0) & (data[tool] != '0') & (data[tool].notna())]
        filtered_tables.append(filtered_data)

    data_filt = pd.concat(filtered_tables).drop_duplicates().reset_index(drop=True)
    
    return data_filt

def amps_in_more_than_x_tools(data, amp_tools, min_tool_n):
    """
    Retains clusters that contain one or more AMP identified by
    X number of tools.
    """
    amp_tools = list(amp_tools)
    
    # Check if 'HMM_model' is in amp_tools and replace it with 'evalue_hmmer'
    if 'HMM_model' in amp_tools:
        amp_tools[amp_tools.index('HMM_model')] = 'HMMsearch'
    
    # Ensure columns are numeric
    for col in amp_tools:
        if col in data.columns:
            data[col] = pd.to_numeric(data[col], errors='coerce')  # Convert to numeric, replace non-numeric with NaN
    # Use replace() in the condition only, without modifying the original DataFrame
    more_than_3 = data[(data[amp_tools] > 0.0).sum(axis=1) >= min_tool_n]
    # Grab the unique clusters that these AMPs belong to
    clusters_3 = more_than_3['cluster_id'].unique()
    # Filter the original DataFrame to include only those clusters
    df_filter = data[data['cluster_id'].isin(clusters_3)]
    
    return df_filter

def amps_with_signalpeptide(data):
    """
    Retain clusters with at least one AMP with a signal peptide
    predicted
    """
    clusters = data[['cluster_id','signal_peptide']]
    # If cluster_id has no 'yes' then remove the dictionary with that cluster
    result = clusters.groupby('cluster_id').agg(signal_peptide=('signal_peptide', list)).reset_index().to_dict('records')
    signal_peptide_list = []
    for d in result:
        if 'yes' in d['signal_peptide']:
            signal_peptide_list.append(d)
    # Retain only AMPs found in filtered clusters
    retain_clusters = [d['cluster_id'] for d in signal_peptide_list]
    df_filter_SP = data[data['cluster_id'].isin(retain_clusters)]
    
    return df_filter_SP