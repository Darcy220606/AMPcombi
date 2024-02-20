import pandas as pd

from typing import Callable
from shiny import Inputs, Outputs, Session, module, render, ui
from shinywidgets import output_widget, render_widget

from plots import plot_ampcombi_upset, plot_ampcombi_cluster, plot_ampcombi_sankey, plot_ampcombi_pdb

####################
#SHINY VERSION == 0.7.1
####################

####################
#AMPCOMBI - TABLE
####################
@module.ui
def ampcombi_table_ui():
    return  ui.nav_panel(
        "Table",
                # download rows selected : table tab
                ui.p(
                    """
                    Download selected rows:
                    """),
                ui.card(
                    ui.download_button("download_ampcombi_table_rows", "Download 'ampcombi_table_selected_rows.tsv'"),
                ),
                # the rows user selected
                {"style": "font-size: 20px;"},
                "Rows selected by user: ", 
                ui.output_text("ampcombi_table_rows", inline=True),
                ui.row(),
                ui.output_data_frame("ampcombi_table_dataframe"),
    )
    
@module.server
def ampcombi_table_server(
    input: Inputs,
    output: Outputs,
    session: Session,
    df: Callable[[], pd.DataFrame],
):
    @render.data_frame
    def ampcombi_table_dataframe():
        """"
        AMPCOMBI: render dataframe in a table
        """
        if isinstance(df(), pd.DataFrame):
            # render grid table
            data_grid = render.DataGrid(df(),                                           
                                            row_selection_mode= 'multiple',width="100%",height="1000px",
                                            filters=True)
            return data_grid
        else: None

    @render.text
    def ampcombi_table_rows():
        """
        AMPCOMBI: prints the row numbers selected by user
        """
        selected = input.ampcombi_table_dataframe_selected_rows()
        l =', '.join(str(i) for i in selected)
        return l

    @render.download(
        filename=lambda: "ampcombi_table_selected_rows.tsv"
    )
    async def download_ampcombi_table_rows():
        indices = list(input.ampcombi_table_dataframe_selected_rows())
        selected_rows = df().iloc[indices]
        yield selected_rows.to_csv('ampcombi_table_selected_rows.tsv', sep='\t', index=False)

    
####################
#AMPCOMBI - UPSET
####################
@module.ui
def ampcombi_upset_ui():
    return  ui.nav_panel(
        "Tool comparison",
                ui.output_plot("ampcombi_upset_plot"),
                #{"style": "font-size: 10px;"},
                #"Note: The upset plot is automatically downloaded as 'ampcombi_upset_plot.pdf/.png'.", 
                ui.output_text("ampcombi_upset_text_download", inline=True),
                ui.p(
                    {"style": "font-size: 20px;"},
                    """
                    Filter the table according to the prediction tool
                    """),
                ui.input_checkbox_group(
                    "ampcombi_upset_checkbox", "", 
                    {"prob_ampir": "Ampir", 
                     "prob_neubi": "Neubi",
                     "HMMsearch":"HMMsearch",
                     "prob_ampgram":"Ampgram",
                     "prob_amptransformer":"Amptransformer",
                     "prob_macrel":"Macrel",
                     "prob_amplify":"Amplify"
                        }
                    ),
                ui.download_button("ampcombi_upset_download_filtered", "Download 'ampcombi_upset_filtered.tsv'"),
                ui.output_data_frame("ampcombi_upset_table"),
    )
    

@module.server
def ampcombi_upset_server(
    input: Inputs,
    output: Outputs,
    session: Session,
    df: Callable[[], pd.DataFrame],
):
    
    @render.plot  
    def ampcombi_upset_plot():
        return plot_ampcombi_upset(df())

    def ampcombi_table_filter_upset():
        if isinstance(df(), pd.DataFrame):
            if input.ampcombi_upset_checkbox() is not None:
                tools_selected = input.ampcombi_upset_checkbox()
                # for shared hits in more than one tool
                if len(tools_selected) == 1 and 'HMMsearch' not in tools_selected:
                    tools = list(tools_selected)
                    filtered_df = df()[(df()[tools] > 0).all(axis=1) & (df()[tools].notna().all(axis=1))]
                    return filtered_df
                elif len(tools_selected) == 1 and 'HMMsearch' in tools_selected:
                    df()['HMM_model'] = df()['HMM_model'].astype('category')
                    filtered_df = df()[(~df()['HMM_model'].isin(['0', '', '0.0'])) & (df()['HMM_model'].notna())]
                    return filtered_df
                elif len(tools_selected) > 1 and 'HMMsearch' not in tools_selected:
                    tools = list(tools_selected)
                    filtered_df = df()[(df()[tools] > 0).all(axis=1) & (df()[tools].notna().all(axis=1))]
                    return filtered_df
                elif len(tools_selected) > 1 and 'HMMsearch' in tools_selected:
                    tools = list(tools_selected)
                    # remove hmmsearch from the list of all tools as that is teh only categorical value in tools
                    filtered_values = list(filter(lambda x: x != 'HMMsearch', tools))
                    filtered_df = df()[(df()[filtered_values] > 0).all(axis=1) & (df()[filtered_values].notna().all(axis=1)) & (~df()['HMM_model'].isin(['0', '', '0.0'])) & (df()['HMM_model'].notna())]
                    return filtered_df
        else:
            None
    
    @render.data_frame
    def ampcombi_upset_table():
        df_amp_filt = ampcombi_table_filter_upset()
        if isinstance(df_amp_filt, pd.DataFrame):
            filtered_table = render.DataTable(df_amp_filt, width="100%")
            return filtered_table
        else:
            return render.DataTable(df(),width="100%")
    
    @render.download(
    filename=lambda: "ampcombi_upset_filtered.tsv"
    )
    async def ampcombi_upset_download_filtered():
        df_amp_filt = ampcombi_table_filter_upset()
        if isinstance(df_amp_filt, pd.DataFrame):
            yield df_amp_filt.to_csv('ampcombi_upset_filtered.tsv', sep='\t', index=False)

####################
#AMPCOMBI - CLUSTER
####################
@module.ui
def ampcombi_cluster_ui():
    return  ui.nav_panel(
        "Clusters",
                output_widget("ampcombi_cluster_plot"),
                ui.p(
                    {"style": "font-size: 20px;"},
                    """
                    Filter the table according to the cluster ID
                    """),
                ui.input_numeric("clusters_id", "", value=1455),
                ui.download_button("ampcombi_cluster_download_filtered", "Download 'ampcombi_cluster_filtered.tsv'"),
                ui.row(),
                ui.output_data_frame("ampcombi_cluster_table"),
    )

@module.server
def ampcombi_cluster_server(
    input: Inputs,
    output: Outputs,
    session: Session,
    df: Callable[[], pd.DataFrame],
):
    
    @render_widget
    def ampcombi_cluster_plot():
        if isinstance(df(), pd.DataFrame):
            return plot_ampcombi_cluster(df())
    
    def ampcombi_table_filter_cluster():
        if isinstance(df(), pd.DataFrame):
            df_cluster = df()[df()['cluster_id'] == input.clusters_id()]
            return df_cluster
        else:
            None

    @render.data_frame
    def ampcombi_cluster_table():
        df_amp_filter = ampcombi_table_filter_cluster()
        if isinstance(df_amp_filter, pd.DataFrame):
            data_grid = render.DataTable(df_amp_filter, width = "100%")
            return data_grid
        else: render.DataTable(df(),width="100%")

    @render.download(
    filename=lambda: "ampcombi_cluster_filtered.tsv"
    )
    async def ampcombi_cluster_download_filtered():
        df_amp_filter = ampcombi_cluster_table()
        if isinstance(df_amp_filter, pd.DataFrame):
            yield df_amp_filter.to_csv('ampcombi_cluster_filtered.tsv', sep='\t', index=False)

####################
#AMPCOMBI - SANKEY
####################
@module.ui
def ampcombi_sankey_ui():
    return  ui.nav_panel(
        "Taxonomy",
                ui.p(
                    {"style": "font-size: 20px;"},
                    """
                    Enter cluster ID here to filter plot
                    """
                    ),
                ui.input_numeric("clusters_id_tax", "", value=None),
                ui.output_plot("ampcombi_sankey_plot"), 
                {"style": "font-size: 15px;"},
                "The sankey plot showing the taxonomic lineage is rendered in a new tab in which the user can interactively adjust.", 
    )

@module.server
def ampcombi_sankey_server(
    input: Inputs,
    output: Outputs,
    session: Session,
    df: Callable[[], pd.DataFrame],
):
    @render.plot  
    def ampcombi_sankey_plot():
        if input.clusters_id_tax() is not None:
            filtered = df()[df()['cluster_id'] == input.clusters_id_tax()]
            return plot_ampcombi_sankey(filtered)
        else: 
            return plot_ampcombi_sankey(df())

####################
#AMPCOMBI - PDB
####################
@module.ui
def ampcombi_pdb_ui():
    return  ui.nav_panel(
        "3D-structure",
                ui.input_file("ampcombi_pdb_files", "Choose a PDB file", accept=[".pdb"], multiple=True),
                {"style": "font-size: 15px;"},
                "Please upload one or more '*.pdb' files here. \n The final structure is saved as 'ampcombi_pdb_structure.html'.", 
                output_widget("ampcombi_pdb_parse"),
    )

@module.server
def ampcombi_pdb_server(
    input: Inputs,
    output: Outputs,
    session: Session,
):
    @render_widget
    def ampcombi_pdb_parse():
        file_infos = input.ampcombi_pdb_files() #its a list of dict, each dict one file 'check: https://shinylive.io/py/examples/#file-upload'
        # checks if file is uploaded by the user
        if not file_infos:
            return None 
        # grab the file paths
        pdb_files = [file_info['datapath'] for file_info in file_infos]
        return plot_ampcombi_pdb(pdb_files)
