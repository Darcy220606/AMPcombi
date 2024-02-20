from modules import ampcombi_table_ui, ampcombi_table_server, ampcombi_upset_ui, ampcombi_upset_server, ampcombi_cluster_ui, ampcombi_cluster_server, ampcombi_sankey_ui, ampcombi_sankey_server, ampcombi_pdb_ui, ampcombi_pdb_server
from shiny import App, Inputs, Outputs, Session, reactive, ui, render

import shinyswatch
from pathlib import Path
import pandas as pd
from shiny.types import FileInfo
from shiny.types import ImgData
import os
import re
from io import StringIO
import asyncio

#################
# UI: user interface function
#################
app_ui = ui.page_navbar(
    shinyswatch.theme.solar(), #darkly/cyborg/superherop
    ampcombi_table_ui("tab1"),
    ampcombi_upset_ui("tab2"),
    ampcombi_cluster_ui("tab3"),
    ampcombi_sankey_ui("tab4"),
    ampcombi_pdb_ui("tab5"),
    sidebar=ui.sidebar(
        # add logo
        #ui.output_image("ampcombi_logo"),
        ui.img(src="amp-combi-logo.png", style="width:200px;"),
            {"style": "color: grey;"},
            {"style": "font-weight: bold;"},
        ui.a(dict(href="https://github.com/Darcy220606/AMPcombi"), "AMPcombi documentation"),
        # upload file in tsv format
        ui.p(
            """
            Choose a file to upload:
            """),
        ui.input_file("ampcombi_user_tsv", label = "",accept=[".tsv"]),     
    width="300px",
    ),
    title="AMPcombi",
    id="tabs",
)

def server(input: Inputs, output: Outputs, session: Session):
    @reactive.Calc()
    def filtered_data() -> pd.DataFrame:
        """
        User uploads table and stores it in function for ampcombi
        """
        file_infos = input.ampcombi_user_tsv() #its a list of dict, each dict one file 'check: https://shinylive.io/py/examples/#file-upload'
        # checks if file is uploaded by the user
        if not file_infos:
            return None 

        out_str = ""
        for file_info in file_infos:
            out_str = file_info["datapath"]
            data = pd.read_csv(out_str, sep='\t')
            #######################################
            # remove unecessary seq_headerss:
            if 'seq_headerss' in data.columns:
                # reorder the columns
                desired_order = list(range(11)) + [69, 70] + list(range(11, 68)) + list(range(71, len(data.columns)))
                data = data.iloc[:, desired_order]
                data['cluster_id'] = data['cluster_id'].astype(int)
            if 'cluster_id' in data.columns:
                data['cluster_id'] = data['cluster_id'].astype(int)
            if 'evalue_hmmer' in data.columns:
                data['evalue_hmmer'] = data['evalue_hmmer'].astype(float)
                data['evalue_hmmer'] = data['evalue_hmmer'].apply(lambda x: '{:.2e}'.format(x))
            if 'HMM_model' in data.columns:
                data['HMM_model'] = data['HMM_model'].astype('category')
            #######################################
        return data
    
    ampcombi_table_server(id="tab1", df=filtered_data)
    ampcombi_upset_server(id="tab2", df=filtered_data)
    ampcombi_cluster_server(id="tab3", df=filtered_data)
    ampcombi_sankey_server(id="tab4", df=filtered_data)
    ampcombi_pdb_server(id="tab5")

# add path to logo
www_dir = Path(__file__).parent / "" #change strings "" to path of where the images should be found ex. "https://github.com/Darcy220606/AMPcombi/blob/main/docs/"
app = App(app_ui, server, static_assets=www_dir)

