import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from upsetplot import UpSet
from typing import Callable
from shiny import Inputs, Outputs, Session, module, render, ui, reactive, req

####################################################
#AMPCOMBI - Upset plot
####################################################
@module.ui
def ampcombi_upset_ui():
    return  ui.nav_panel(
        "Upset",
        ui.card(
                ui.output_plot("ampcombi_upset_plot")),
        ui.card(
                ui.p("Table:"),
                ui.output_data_frame("ampcombi_table_upset"),
                )
    )

@module.server
def ampcombi_upset_server(
    input: Inputs,
    output: Outputs,
    session: Session,
    df: Callable[[], pd.DataFrame]
    ):
    
    @output
    @render.data_frame
    def ampcombi_table_upset():
        """"
        AMPCOMBI: render dataframe in a table
        """
        if isinstance(df(), pd.DataFrame):
            # Render grid table
            return render.DataTable(df(), width="100%")
        else: None    

    @output
    #@render.ui
    @render.plot
    def ampcombi_upset_plot():
        if isinstance(df(), pd.DataFrame):
            return upset_plot(df())
        else: None
    
# Create upset plot from data
def upset_plot(table):
    """
    Adjusts the input for the upset plot.
    
    """
    # Define the tools
    #column_order=["prob_ampir","prob_ampgram","prob_neubi","prob_amptransformer","HMM_model","prob_amplify","prob_macrel"]
    column_order=["AMPir","AMPgram","NEUBI","AMPtransformer","HMM_model","AMPlify","MACREL"]
    # Check if any values from the list are present as columns
    matching_columns = [col for col in column_order if col in table.columns]
    # If matching columns exist, create a new DataFrame
    if matching_columns:
        table_df = table[matching_columns]
    # Replace 0.0 with NAs
    table_df.replace(0.0, np.nan, inplace=True)
    table_df.replace(0, np.nan, inplace=True)
    # Turn df into 0s and 1s
    table_df = table_df.applymap(lambda x: 0 if pd.isna(x) else 1)
    # Convert to true false 
    table_df = table_df.astype(bool)
    # Reformat for upset plot
    multi_index_series = table_df.groupby(list(table_df.columns)).size()

    # Recolour
    tool_color_dict = {
        'AMPir': '#E6AB02',     # Mustard yellow
        'AMPgram': '#A6761D',   # Burnt sienna
        'NEUBI': '#666666',     # Charcoal grey
        'AMPtransformer': '#1B9E77', # Teal
        'MACREL': '#D95F02',    # Orange-red
        'AMPlify': '#7570B3',   # Lavender
        'HMM_model': '#E7298A'       # Fuchsia
    }

    # Only keep colors for the columns that are present
    tool_color_dict_filtered = {key: tool_color_dict[key] for key in matching_columns if key in tool_color_dict}

    # Plot
    fig = plt.figure(figsize=(14, 4), facecolor="white")

    # With plt.style.context('dark_background'):
    upset = UpSet(
        multi_index_series, 
        sort_by='degree', 
        show_counts=True, 
        facecolor="black", #lightblue
        other_dots_color=.3, 
        shading_color=.1, 
        totals_plot_elements=4, 
        element_size=20) # reduce dot size
    upset_plot = upset.plot()

    # Set axes for plot on bottom left corner with totals of each set-
    set_bar_ax = upset_plot["totals"]
    # Adjust the colour to each set 
    for patch, label in zip(set_bar_ax.patches, table_df.columns):
        if label in tool_color_dict:
            patch.set_facecolor(tool_color_dict_filtered[label])
    
    #fig.savefig("debug_output.png", format="png", bbox_inches="tight")
    