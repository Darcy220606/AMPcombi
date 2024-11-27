import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px

from shinywidgets import output_widget, render_widget
from typing import Callable
from shiny import Inputs, Outputs, Session, module, render, ui, reactive, req

####################################################
#AMPCOMBI - HMM models
####################################################

@module.ui
def ampcombi_ampclasses_ui():
    return  ui.nav_panel(
        "AMP-class",
        ui.input_select("ref_db_amp_class", "Select the database used to classify AMP classes:", 
                            {"HMM_model":"HMM models", 
                             "interpro_accession": "InterProScan",
                             "DRAMP_ID": "DRAMP",
                             "APD_ID" : "APD",
                             "header":"UniRef100"}, selected=None),
        ui.card(
            ui.p("AMP classes across the tools:"),
            output_widget("refDB_stacked_tool")),
        ui.card(
            ui.p("AMP classes across the samples:"),
            output_widget("refDB_stacked_sample")),
    )

@module.server
def ampcombi_ampclasses_server(
    input: Inputs,
    output: Outputs,
    session: Session,
    df: Callable[[], pd.DataFrame]
    ):
    
    @render_widget
    def refDB_stacked_tool():
        """
        Plot of the proportions of AMP classes per tool.
        """
        if isinstance(df(), pd.DataFrame):
            return refDB_stacked_plot_tool(df(), input.ref_db_amp_class())

    @render_widget
    def refDB_stacked_sample():
        """
        Plot of the proportions of AMP classes per sample.
        """
        if isinstance(df(), pd.DataFrame):
            return refDB_stacked_plot_sample(df(), input.ref_db_amp_class())

def refDB_stacked_plot_tool(data, model_column):
    """
    Create a stacked bar plot for AMP prediction tools and a specified model column.
    """
    if model_column not in data.columns:
        raise ValueError(f"The specified column '{model_column}' does not exist in the dataset. "
                         f"Please use an alternative database.")
    # Replace 0.0 with NAs
    data.replace(0.0, np.nan, inplace=True)
            
    def process_tool(tool_name, data, model_column):
        """
        Process data for a specific tool and model column.
        """
        tool_df = data[[tool_name, model_column]]
        tool_df = tool_df.dropna(subset=[tool_name])
        tool_counts = tool_df[model_column].value_counts()
        tool_df_2 = pd.DataFrame({model_column: tool_counts.index, tool_name: tool_counts.values})
        #tool_df_2.iloc[0,0] == 'Unknown'
        return tool_df_2
    
    # Identify available tools dynamically
    tool_columns = ["AMPir", "AMPgram", "NEUBI", "AMPtransformer", "AMPlify", "MACREL", "HMMsearch"]
    available_tools = [tool for tool in tool_columns if tool in data.columns]
    
    tool_dfs = []
    for tool in available_tools:
        
        tool_dfs.append(process_tool(tool, data, model_column))

    # Merge all processed DataFrames horizontally
    if tool_dfs:
        final_df = tool_dfs[0]
        for tool_df in tool_dfs[1:]:
            final_df = final_df.merge(tool_df, on=model_column, how="outer")

    # Convert to long format
    final_df_long = pd.melt(final_df, id_vars=[model_column], var_name="Tool", value_name="Counts")
    
    # Plot using Plotly
    fig = px.bar(
        final_df_long,
        x="Tool",
        y="Counts",
        color=model_column,
        #title=f"AMP Prediction Tools vs {model_column}",
        color_discrete_sequence=px.colors.sequential.Viridis_r,
        labels={"Tool": "AMP Prediction Tools", "Counts": "Number of AMPs"},
    )

    # Customize layout
    fig.update_layout(
        xaxis_title="AMP Prediction Tools",
        yaxis_title="Number of AMPs",
        legend_title=model_column,
        barmode="stack",
        title_text="Database IDs",
        #width=800,
        #height=500,
        legend=dict(
            orientation="v",
            y=1,
            x=1.02,
            xanchor="left",
            yanchor="top",
        ),
    )

    # Rotate x-axis labels for better visibility
    return fig.update_xaxes(tickangle=25)

def refDB_stacked_plot_sample(data, model_column):
    """
    Create a stacked bar plot for AMP prediction tools and a specified model column.
    """
    if model_column not in data.columns:
        raise ValueError(f"The specified column '{model_column}' does not exist in the dataset. "
                         f"Please use an alternative database.")
    # Replace 0.0 with NAs
    data.replace(0.0, np.nan, inplace=True)
    
    # Count AMPs per sample and per model
    counts = data.groupby(["sample_id", model_column]).size().reset_index(name="Counts")

    # Plot using Plotly
    fig = px.bar(
        counts,
        x="sample_id",
        y="Counts",
        color=model_column,
        title="AMP Counts Per Sample by Model",
        color_discrete_sequence=px.colors.sequential.Viridis_r,
        labels={"sample_id": "Sample ID", "Counts": "Number of AMPs"},
    )

    # Customize layout
    fig.update_layout(
        xaxis_title="Sample ID",
        yaxis_title="Number of AMPs",
        legend_title=model_column,
        barmode="stack",
        #title_text="AMP Counts Per Sample by Model",
        legend=dict(
            orientation="v",
            y=1,
            x=1.02,
            xanchor="left",
            yanchor="top",
        ),
    )

    # Rotate x-axis labels for better visibility
    return fig.update_xaxes(tickangle=25)

