import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px

from shinywidgets import output_widget, render_widget
from typing import Callable
from shiny import Inputs, Outputs, Session, module, render, ui, reactive, req

####################################################
#AMPCOMBI - Stats
####################################################
@module.ui
def ampcombi_stats_ui():
    return  ui.nav_panel(
        "Stats",
        ui.markdown("All graphs pertaining to the physio-chemical properties."),
        ui.card(
            ui.p("AMP Prediction Score vs Length by Tool:"),
            output_widget("lengths")),
            ui.card(
            ui.p("AMP Properties by Tool:"),
            ui.input_select(
                "xaxis_stats",
                "Choose x-axis:",
                {"AMP Length":"AMP lengths",
                 "MwT": "MwT","IP": "IP", "Hydrophobicity" : "Hydrophobicity",
                 "Helix fraction": "Helix proportion",  "Sheet fraction": "Sheet proportion",  
                 "Turn fraction": "Turn proportion"},),
            output_widget("physio_chemical_props"),
            ),
    )

@module.server
def ampcombi_stats_server(
    input: Inputs,
    output: Outputs,
    session: Session,
    df: Callable[[], pd.DataFrame]
    ):
    
    #@output
    @render_widget
    def lengths():
        """"
        Plot of the distribution of the lengths of the molecules.
        """
        if isinstance(df(), pd.DataFrame):
            return len_vs_score(df())
        else: None    
        
    @render_widget
    def physio_chemical_props():
        """
        Plot of the distribution of the physico-chemical properties of the prepropeptides.
        """
        if isinstance(df(), pd.DataFrame):
            return x_vs_y_axis_counts(df(), input.xaxis_stats())

# Create regression plot from data
def len_vs_score(data):
    """
    Plot scatter plots with lengths VS score
    """
    # Define the tools
    column_order=["CDS_id", "AMP Length", "AMPir","AMPgram","NEUBI","AMPtransformer","AMPlify","MACREL"]
    # Check if any values from the list are present as columns
    matching_columns = [col for col in column_order if col in data.columns]
    # If matching columns exist, create a new DataFrame
    if matching_columns:
        table_df = data[matching_columns]

    unfiltered_long = pd.melt(table_df, id_vars=['CDS_id','AMP Length'], var_name='Tool',value_name='Score')
    unfiltered_long['Score'] = unfiltered_long['Score'].replace(0, np.nan)
    
    tool_color_dict = {
        'AMPir': '#E6AB02',     # Mustard yellow
        'AMPgram': '#A6761D',   # Burnt sienna
        'NEUBI': '#666666',     # Charcoal grey
        'AMPtransformer': '#1B9E77', # Teal
        'MACREL': '#D95F02',    # Orange-red
        'AMPlify': '#7570B3',   # Lavender
        'HMM_model': '#E7298A'       # Fuchsia
    }

    # Create the scatter plot with facets
    fig = px.scatter(
        unfiltered_long,
        x="AMP Length", 
        y="Score",      
        color="Tool",   
        color_discrete_map=tool_color_dict,
        facet_col="Tool",
        facet_col_wrap=3,
        trendline="ols", 
        #title="AMP Prediction Score vs Length by Tool"
    )
    
    return fig.update_layout(
    xaxis_title="AMP Length",
    yaxis_title="Prediction Score",
    title_font_size=16,
    font=dict(size=12),
    height=600,
    width=1000)

# Create regression plot from data
def x_vs_y_axis_counts(data, xaxis):
    """
    Plot scatter plots with X (diff physiochemical properties) against the number of AMPS.
    """     
    # Hydrophobicity special case
    data['Hydrophobicity'].replace('None', 0.0, inplace=True)
    data['Hydrophobicity'] = data['Hydrophobicity'].astype(float)
    # Replace 0s with NAs
    data.replace(0.0, np.nan, inplace=True)
    
    # Identify available tools dynamically
    tool_columns = ["AMPir", "AMPgram", "NEUBI", "AMPtransformer", "AMPlify", "MACREL", "HMMsearch"]
    available_tools = [tool for tool in tool_columns if tool in data.columns]

    # Function to create the DataFrame for each tool
    def process_tool_xaxis(tool_name, df, xaxis):
        tool_df = df[[tool_name, xaxis]]
        tool_df = tool_df.dropna(subset=[tool_name])
        tool_counts = tool_df[xaxis].value_counts()
        tool_df_2 = pd.DataFrame({xaxis: tool_counts.index, tool_name: tool_counts.values})
        tool_df_2 = tool_df_2.sort_values(tool_name, ascending=False)
        return tool_df_2

    # Process only available tools
    tool_dfs = []
    for tool in available_tools:
        tool_dfs.append(process_tool_xaxis(tool, data, xaxis))

    # Merge all processed DataFrames horizontally
    if tool_dfs:
        final_df = tool_dfs[0]
        for tool_df in tool_dfs[1:]:
            final_df = final_df.merge(tool_df, on=xaxis, how="outer")

    # Convert to long format
    final_df_long = pd.melt(final_df, id_vars=[xaxis], var_name="Tool", value_name="Counts")
    
    tool_color_dict = {
        'AMPir': '#E6AB02',     # Mustard yellow
        'AMPgram': '#A6761D',   # Burnt sienna
        'NEUBI': '#666666',     # Charcoal grey
        'AMPtransformer': '#1B9E77', # Teal
        'MACREL': '#D95F02',    # Orange-red
        'AMPlify': '#7570B3',   # Lavender
        'HMM_model': '#E7298A'       # Fuchsia
    }

    # Plot the data
    fig = px.scatter(
        final_df_long,
        x=xaxis, y="Counts",
        trendline="lowess",
        color="Tool",   
        color_discrete_map=tool_color_dict,
        facet_col="Tool",
        facet_col_wrap=3,
    )

    # Update layout
    return fig.update_layout(
        xaxis_title=xaxis,
        yaxis_title="Number of AMPs",
        title_font_size=16,
        font=dict(size=12),
        height=600,
        width=1000,
    )
