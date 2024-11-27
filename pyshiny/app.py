
# MIT License
# 
# Copyright (c) 2024 Anan Ibrahim (darcy220606)
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from grid import (ampcombi_table_ui, ampcombi_table_server, datatable_filter)
from upset import (ampcombi_upset_ui, ampcombi_upset_server)
from stats import (ampcombi_stats_ui, ampcombi_stats_server)
from ampclasses import (ampcombi_ampclasses_ui, ampcombi_ampclasses_server)
from clusters import (ampcombi_cluster_ui, ampcombi_cluster_server)
from sankey import (ampcombi_sankey_ui, ampcombi_sankey_server)

from shiny import App, Inputs, Outputs, Session, reactive, ui, render, req

from pathlib import Path
import pandas as pd
from shiny.types import FileInfo
from shiny.types import ImgData
import os
import re

#################
# UI: user interface function
#################
# Main page structure
app_ui = ui.page_navbar(
    ui.nav_spacer(), 
    ui.nav_control(ui.input_dark_mode()), 
    ampcombi_table_ui("tab1"),
    ampcombi_upset_ui("tab2"),
    ampcombi_stats_ui("tab3"),
    ampcombi_ampclasses_ui("tab4"),
    ampcombi_cluster_ui("tab5"),
    ampcombi_sankey_ui("tab6"),
    sidebar=ui.sidebar(
        # Add logo
        ui.img(src="amp-combi-logo.png", style="width:200px;"),
        ui.a(dict(href="https://ampcombi.readthedocs.io/en/main/"), "AMPcombi documentation"),
        # Upload file in tsv format
        ui.p("Choose a file to upload:"),
        ui.input_file(
            "ampcombi_user_tsv", 
            label = "",
            accept=[".tsv"]),
        ui.HTML("<h4 style='color: #808080; font-size: 18px; font-weight: bold; margin-bottom: -5px;'>Filter by AMPs</h4>"),
        # Checkbox tools
        ui.p("Select Prediction Tool:"),
        #ui.HTML("<h4 style='color: #595959; font-size: 18px; font-weight: regular; margin-bottom: -5px;'>Select Prediction Tool</h4>"),
        ui.input_checkbox_group(
                "tool_selection", 
                None,
                {"AMPir": "ampir", 
                 "NEUBI": "neubi",
                 "AMPgram":"ampgram",
                 "MACREL":"macrel",
                 "AMPlify":"amplify",
                 "AMPtransformer":"amptransformer",
                 "HMM_model":"HMM model"
                },
                selected=["AMPir","AMPgram","NEUBI","AMPtransformer","HMM_model","AMPlify","MACREL"]),
                #choices=["ampir", "amplify", "ampgram", "neubi", "macrel", "amptransformer", "HMMsearch"], 
        ui.input_slider("slider_prediction_score", "Min. Prediction score", min=0, max=1, value=0), 
        ui.input_slider("slider_amp_lengths", "Max. AMP lengths",  min=1, max=2000, value=1500), 
        ui.HTML("<h4 style='color: #808080; font-size: 18px; font-weight: bold; margin-bottom: -5px;'>Filter by Clusters</h4>"),
        ui.input_selectize("pick_clusters", "Choose a cluster:", [], multiple=True),
        ui.row(
            ui.p("Retain clusters with at least X AMP identified with:",style="font-size: 15px;"),
            ui.input_numeric("min_tool_number", "Min. tool number", 1, min=1, max=5),
            ui.input_switch("signal_peptide", "Remove clusters with no signaling peptides", False),
            ), 
    ),   
    title="AMPcombi",
    id="tabs",
)

#################
# Backend: server + reactivity
#################
def server(input: Inputs, output: Outputs, session: Session):
    @reactive.Calc()
    def filtered_data() -> pd.DataFrame:
        """
        User uploads table and stores it in function for ampcombi
        """
        file_infos = input.ampcombi_user_tsv() #its a list of dict, each dict one file 'check: https://shinylive.io/py/examples/#file-upload'
        # Checks if file is uploaded by the user
        if not file_infos:
            return None 

        for file_info in file_infos:
            out_str = file_info["datapath"]
            data = pd.read_csv(out_str, sep='\t')
        
        # Grab the side selections
        amp_tools = input.tool_selection()
        prediction_score = input.slider_prediction_score()
        amp_lengths = input.slider_amp_lengths()
        min_tool_n = int(input.min_tool_number())
        remove_signal_peptide = input.signal_peptide()
        
        # Apply the filter_data function to filter based on sidebar inputs
        data = datatable_filter(data, amp_tools, prediction_score, amp_lengths, min_tool_n, remove_signal_peptide)

        # Extract unique cluster IDs if `cluster_id` exists
        if "cluster_id" in data.columns:
            unique_clusters = sorted(data["cluster_id"].unique().tolist())  # Sort alphabetically
            
            # Preserve current selection
            selected_clusters = input.pick_clusters()
            # Dynamically update selectize options while preserving the current selection
            ui.update_selectize(
                "pick_clusters",
                choices=unique_clusters,
                selected=selected_clusters,
                server=False,  # Enable server-side updates for performance
            )
            if selected_clusters:
                data = data[data["cluster_id"].isin(selected_clusters)]
        
        return data
    
    ampcombi_table_server(id="tab1", df=filtered_data)
    ampcombi_upset_server(id="tab2", df=filtered_data)
    ampcombi_stats_server(id="tab3", df=filtered_data)
    ampcombi_ampclasses_server(id="tab4", df=filtered_data)
    ampcombi_cluster_server(id="tab5", df=filtered_data)
    ampcombi_sankey_server(id="tab6", df=filtered_data)

# Add path to logo
www_dir = Path(__file__).parent / "" #change strings "" to path of where the images should be found ex. "https://github.com/Darcy220606/AMPcombi/blob/main/docs/"
app = App(app_ui, server, static_assets=www_dir)

