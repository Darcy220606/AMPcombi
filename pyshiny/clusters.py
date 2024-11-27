import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

from shinywidgets import output_widget, render_widget
from typing import Callable
from shiny import Inputs, Outputs, Session, module, render, ui, reactive, req
from sklearn.decomposition import PCA
#from sklearn.metrics.pairwise import cosine_similarity

####################################################
#AMPCOMBI - HMM models
####################################################
@module.ui
def ampcombi_cluster_ui():
    return  ui.nav_panel(
        "AMP-clusters",
        ui.input_select("ref_db_amp_class", "Select the database used to classify AMP classes:", 
                            {"HMM_model":"HMM models", 
                             "interpro_accession": "InterProScan",
                             "DRAMP_ID": "DRAMP",
                             "APD_ID" : "APD",
                             "header":"UniRef100"}, selected=None),
        ui.card(
            ui.p("AMP clusters across the samples:"),
            output_widget("amp_clusters_samples")),
        ui.card(
            ui.p("AMP clusters across taxonomy:"),
            ui.input_selectize(
                "pick_taxa_group",
                "Choose a taxonomic level:",
                choices=['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'specie'],
                selected='kingdom'
            ),
            output_widget("amp_clusters_taxa")),
        ui.card(
            ui.p("Sample similarity based on gene clusters:"),
            output_widget("sample_similarity")),
    )

@module.server
def ampcombi_cluster_server(
    input: Inputs,
    output: Outputs,
    session: Session,
    df: Callable[[], pd.DataFrame]
    ):
    
    @render_widget
    def amp_clusters_samples():
        """
        AMP clusters across the samples
        """
        if isinstance(df(), pd.DataFrame):
            return amp_clusters_per_sample(df(), input.ref_db_amp_class())

    @render_widget
    def amp_clusters_taxa():
        """
        AMP clusters across the taxa
        """
        if isinstance(df(), pd.DataFrame):
            return amp_clusters_per_taxa(df(), input.ref_db_amp_class(), input.pick_taxa_group())

    @render_widget
    def sample_similarity():
        """
        Sample similarity using teh cluster AMP content.
        """
        if isinstance(df(), pd.DataFrame):
            return sample_similarity_with_clusters(df())

def amp_clusters_per_sample(data, model_column):
    """
    Generate a heatmap of AMP clusters per sample using Plotly.
    """
    if model_column not in data.columns:
        return None    
    if 'cluster_id' not in data.columns:
        return None

    # Replace 0.0 with NaNs
    data.replace(0.0, np.nan, inplace=True)

    # Set up the input data
    sample_df = data[['sample_id', 'cluster_id']]

    # Long to short format
    sample_df = sample_df.groupby(['sample_id', 'cluster_id'], dropna=False).size().reset_index(name='count')
    df_wide = sample_df.pivot(index='cluster_id', columns='sample_id', values='count')

    # Sort by presence across all samples
    df_wide['abund'] = df_wide.sum(axis=1)
    df_sorted = df_wide.sort_values(by=['abund'], ascending=[False])
    df_sorted = df_sorted.drop(columns=['abund'])
    df_sorted.fillna(0, inplace=True)  # Replace NaNs with 0
    df_sorted = df_sorted.astype(int)

    if model_column == 'DRAMP_ID':
        refDB_sp_col = 'Name'
    elif model_column == 'APD_ID':
        refDB_sp_col = 'target'
    elif  model_column == 'HMM_model':
        refDB_sp_col = 'HMM_model'
    elif  model_column == 'interpro_accession':
        refDB_sp_col = 'interpro_description'
    elif model_column == 'header':
        refDB_sp_col = 'header'

    # Replace 0 with NaN
    data[refDB_sp_col] = data[refDB_sp_col].replace(0, np.nan)
    data[refDB_sp_col] = data[refDB_sp_col].replace('0', np.nan)

    # Set up a dictionary with cluster ID and the list of possible AMP classes
    label_dict = data.groupby('cluster_id')[refDB_sp_col].apply(lambda x: x[pd.notna(x)].tolist()).to_dict()  # Group and ignore NaN
    #label_dict = {key: list(set(value)) for key, value in label_dict.items()}  # Remove duplicates in the list in the dictionary
    label_dict = {key: sorted(set(value)) for key, value in label_dict.items()}

    # Prepare data
    z = df_sorted.values 
    x = df_sorted.columns.tolist()
    y = df_sorted.index.tolist()  

    # Add cluster AMP class labels as annotations
    annotations = []
    max_annotation_length = 0  # Track the longest annotation for margin adjustment

    for i, cluster_id in enumerate(y):
        if cluster_id in label_dict:
            # Join AMP classes into a single string
            label_text = ', '.join(label_dict[cluster_id])  # Combine AMP classes
            max_annotation_length = max(max_annotation_length, len(label_text))

            # Create the annotation and append to the list
            annotations.append(
                dict(
                    x=len(x) + 0.5,  # Position outside the heatmap
                    y=i,  # Align with the cluster ID row
                    text=label_text,  # Keep it as a single string
                    showarrow=False,
                    font=dict(size=10, color="black"),
                    xanchor="left",
                    align="left",  # Align text to the left
                )
            )

    # Prepare the hover text
    hover_text = []
    for i, cluster_id in enumerate(y):
        row_text = []
        for j, sample_id in enumerate(x):
            # Check if the cell value is non-NaN
            if not np.isnan(z[i][j]):
                if cluster_id in label_dict:
                    # AMP class names, split onto multiple lines
                    label_text = '<br>'.join(label_dict[cluster_id])  # Separate classes with HTML line breaks
                    hover_info = (
                        f"Cluster: {cluster_id}<br>"
                        f"Sample: {sample_id}<br>"
                        f"AMPs: {z[i][j]}<br>"
                        f"Annotation:<br>{label_text}"
                    )
                else:
                    hover_info = f"Cluster: {cluster_id}<br>Sample: {sample_id}<br>AMPs: {z[i][j]}"
            else:
                hover_info = ""
            row_text.append(hover_info)
        hover_text.append(row_text)
    
    # Create the Plotly heatmap
    fig = go.Figure()

    # Remove 0 AMPs
    z = np.where(z == 0, None, z)
    
    # Add the heatmap
    fig.add_trace(
        go.Heatmap(
            z=z,
            x=x,
            y=y,
            colorscale="Viridis_r",
            colorbar=dict(title="Number of AMPs", thickness=12),
            hoverinfo="text",  # Use custom hover text
            text=hover_text,   # Assign the hover text
        )
    )

    # Add AMP class labels as annotations
    fig.update_layout(annotations=annotations)

    # Dynamically adjust the right margin based on annotation length
    right_margin = min(max_annotation_length * 6, 500)  # Adjust scale or cap the maximum margin

    # Customize axis titles and layout
    return fig.update_layout(
        xaxis=dict(title="Sample ID", tickangle=45),
        yaxis=dict(
            title="Cluster ID",
            tickmode="array",
            tickvals=list(range(len(y))),  # Ensure all rows are represented
            ticktext=y,  # Use cluster IDs for labels
            tickfont=dict(size=12),  # Increase font size for better visibility
            dtick=1  # Explicitly space rows equally
        ),
        margin=dict(l=150, r=right_margin, t=50, b=50),  # Increase right margin dynamically
        width=1400,  # Set figure width
        height=1000   # Set figure height
    )

def amp_clusters_per_taxa(data, model_column, taxa):
    """
    Plot the number of AMP clusters per taxa.
    """
    # Check input presence
    if 'mmseqs_lineage_contig' in data.columns:
        data["lineage"] = data["mmseqs_lineage_contig"].str.extract(r"(d_.*)")
        data[['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'specie']] = data['mmseqs_lineage_contig'].str.split(';', expand=True)
        # Remove the prefix from each column
        data['kingdom'] = data['kingdom'].str.replace('d_', '')
        data['phylum'] = data['phylum'].str.replace('p_', '')
        data['class'] = data['class'].str.replace('c_', '')
        data['order'] = data['order'].str.replace('o_', '')
        data['family'] = data['family'].str.replace('f_', '')
        data['genus'] = data['genus'].str.replace('g_', '')
        data['specie'] = data['specie'].str.replace('s_', '')
        # Remove the letters used in GTDB formating
        data['phylum'] = data['phylum'].str.replace(r'\s[A-Z](?!\w)', '', regex=True)
        data['class'] = data['class'].str.replace(r'\s[A-Z](?!\w)', '', regex=True)
        data['order'] = data['order'].str.replace(r'\s[A-Z](?!\w)', '', regex=True)
        data['family'] = data['family'].str.replace(r'\s[A-Z](?!\w)', '', regex=True)
        data['genus'] = data['genus'].str.replace(r'\s[A-Z](?!\w)', '', regex=True)
        data['specie'] = data['specie'].str.replace(r'\s[A-Z](?!\w)', '', regex=True)
        # Remove the genus from specie column
        # data['specie_only'] = data['specie'].str.split(' ', n=1).str[1]
        # Replace empty strings with unclassified
        columns_to_update = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'specie']
        data[columns_to_update] = data[columns_to_update].fillna('unclassified')
    else:
        return None
    if model_column not in data.columns:
        return None    
    if 'cluster_id' not in data.columns:
        return None

    # Replace 0.0 with NaNs
    data.replace(0.0, np.nan, inplace=True)

    # Set up the input data
    sample_df = data[[taxa, 'cluster_id']]

    # Long to short format
    sample_df = sample_df.groupby([taxa, 'cluster_id'], dropna=False).size().reset_index(name='count')
    df_wide = sample_df.pivot(index='cluster_id', columns=taxa, values='count')
    # Sort by presence across all samples
    df_wide['abund'] = df_wide.sum(axis=1)
    df_sorted = df_wide.sort_values(by=['abund'], ascending=[False])
    df_sorted = df_sorted.drop(columns=['abund'])
    df_sorted.fillna(0, inplace=True)  # Replace NaNs with 0
    df_sorted = df_sorted.astype(int)

    if model_column == 'DRAMP_ID':
        refDB_sp_col = 'Name'
    elif model_column == 'APD_ID':
        refDB_sp_col = 'target'
    elif  model_column == 'HMM_model':
        refDB_sp_col = 'HMM_model'
    elif  model_column == 'interpro_accession':
        refDB_sp_col = 'interpro_description'
    elif model_column == 'header':
        refDB_sp_col = 'header'

    # Replace 0 with NaN
    data[refDB_sp_col] = data[refDB_sp_col].replace(0, np.nan)
    data[refDB_sp_col] = data[refDB_sp_col].replace('0', np.nan)

    # Set up a dictionary with cluster ID and the list of possible AMP classes
    label_dict = data.groupby('cluster_id')[refDB_sp_col].apply(lambda x: x[pd.notna(x)].tolist()).to_dict()  # Group and ignore NaN
    label_dict = {key: list(set(value)) for key, value in label_dict.items()}  # Remove duplicates in the list in the dictionary

    # Prepare data
    z = df_sorted.values 
    x = df_sorted.columns.tolist()
    y = df_sorted.index.tolist()  

    # Add cluster AMP class labels as annotations
    annotations = []
    max_annotation_length = 0  # Track the longest annotation for margin adjustment

    for i, cluster_id in enumerate(y):
        if cluster_id in label_dict:
            # Join AMP classes into a single string
            label_text = ', '.join(label_dict[cluster_id])  # Combine AMP classes
            max_annotation_length = max(max_annotation_length, len(label_text))

            # Create the annotation and append to the list
            annotations.append(
                dict(
                    x=len(x) + 0.5,  # Position outside the heatmap
                    y=i,  # Align with the cluster ID row
                    text=label_text,  # Keep it as a single string
                    showarrow=False,
                    font=dict(size=10, color="black"),
                    xanchor="left",
                    align="left",  # Align text to the left
                )
            )

    # Prepare the hover text
    hover_text = []
    for i, cluster_id in enumerate(y):
        row_text = []
        for j, sample_id in enumerate(x):
            # Check if the cell value is non-NaN
            if not np.isnan(z[i][j]):
                if cluster_id in label_dict:
                    # AMP class names, split onto multiple lines
                    label_text = '<br>'.join(label_dict[cluster_id])  # Separate classes with HTML line breaks
                    hover_info = (
                        f"Cluster: {cluster_id}<br>"
                        f"Sample: {sample_id}<br>"
                        f"AMPs: {z[i][j]}<br>"
                        f"Annotation:<br>{label_text}"
                    )
                else:
                    hover_info = f"Cluster: {cluster_id}<br>Taxonomy: {sample_id}<br>AMPs: {z[i][j]}"
            else:
                hover_info = ""
            row_text.append(hover_info)
        hover_text.append(row_text)
    
    # Create the Plotly heatmap
    fig = go.Figure()
    
    # Remove 0 AMPs
    z = np.where(z == 0, None, z)

    # Add the heatmap
    fig.add_trace(
        go.Heatmap(
            z=z,
            x=x,
            y=y,
            colorscale="Viridis_r",
            colorbar=dict(title="Number of AMPs", thickness=12),
            hoverinfo="text",  # Use custom hover text
            text=hover_text,   # Assign the hover text
        )
    )

    # Add AMP class labels as annotations
    fig.update_layout(annotations=annotations)

    # Dynamically adjust the right margin based on annotation length
    right_margin = min(max_annotation_length * 6, 500)  #

    # Customize axis titles and layout
    return fig.update_layout(
        xaxis=dict(title="Taxanomic level", tickangle=45),
        yaxis=dict(
            title="Cluster ID",
            tickmode="array",
            tickvals=list(range(len(y))),  # Ensure all rows are represented
            ticktext=y,  # Use cluster IDs for labels
            tickfont=dict(size=12),  # Increase font size for better visibility
            dtick=1  # Explicitly space rows equally
        ),
        margin=dict(l=150, r=right_margin, t=50, b=50),  # Increase right margin dynamically
        width=1400,  # Set figure width
        height=1000   # Set figure height
    )

def sample_similarity_with_clusters(data):
    if 'cluster_id' not in data.columns:
        return None
    # Replace 0.0 with NaNs
    data.replace(0.0, np.nan, inplace=True)

    # Set up the input data
    sample_df = data[['sample_id', 'cluster_id']]

    # Long to short format
    sample_df = sample_df.groupby(['sample_id', 'cluster_id'], dropna=False).size().reset_index(name='count')
    df_wide = sample_df.pivot(index='sample_id', columns='cluster_id', values='count')
    df_wide = df_wide.fillna(0)

    ## Compute similarity matrix
    #similarity_matrix = cosine_similarity(df_wide)
    ## Convert to DataFrame for readability
    #similarity_df = pd.DataFrame(similarity_matrix, index=df_wide.index, columns=df_wide.index)
    #
    ## Create the Plotly heatmap
    #fig = go.Figure(data=go.Heatmap(
    #    z=similarity_df.values,   # Similarity values
    #    x=similarity_df.columns,
    #    y=similarity_df.index,
    #    colorscale="viridis_r",
    #    colorbar=dict(title="Similarity"),  # Colorbar label
    #    text=similarity_df.values.round(2),  # Hover text rounded to 2 decimals
    #    hovertemplate="Similarity between %{x} and %{y}: %{text}<extra></extra>"  # Custom hover template
    #))
#
    ## Customize layout
    #return fig.update_layout(
    #    xaxis=dict(title="Sample ID", tickangle=45),
    #    yaxis=dict(title="Sample ID"),
    #    width=800,  # Adjust width
    #    height=700  # Adjust height
    #)

    # Apply PCA
    pca = PCA(n_components=2)  # Reduce to 2 components
    pca_results = pca.fit_transform(df_wide)

    # Create a DataFrame for PCA results
    pca_df = pd.DataFrame(
        pca_results, 
        columns=["PC1", "PC2"], 
        index=df_wide.index
    ).reset_index()

    # Plot with Plotly
    fig = px.scatter(
        pca_df,
        x="PC1",
        y="PC2",
        text="sample_id",
        labels={"PC1": "Principal Component 1", "PC2": "Principal Component 2"},
        hover_name="sample_id",
        width=800,
        height=600,
    )
    
    return fig.update_traces(marker=dict(size=10, color="blue"), textposition="top center")


        