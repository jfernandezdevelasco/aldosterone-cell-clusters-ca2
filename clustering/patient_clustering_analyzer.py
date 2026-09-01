from typing import Dict, List, Tuple, Union, Optional
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform, pdist

########## Module for analyzing patient single cell calcium trace clustering data.

class PatientAnalyzer:

    @staticmethod
    def calculate_dissimilarity(data):
        # Compute dissimilarity matrix based on Pearson correlation.
        corr_matrix = data.corr()
        return 1.0 - corr_matrix

    #Perform hierarchical agglomerative clustering at a specified distance threshold.
    def predict_clusters(self, data, distance_threshold, method = "complete"): # Complete linkage
        dissimilarity = self.calculate_dissimilarity(data)
        condensed_dist = squareform(dissimilarity.to_numpy(), checks=False)
        linkage_matrix = linkage(condensed_dist, method=method)
        return fcluster(linkage_matrix, distance_threshold, criterion="distance") # prediction at certain 'distance'

    # Extract functional cluster metrics and compare predicted vs true labels.
    def analyze_clusters(self, predicted_labels, true_labels):
        
        df = pd.DataFrame({"Pred": predicted_labels, "Labels": true_labels})
        crosstab = pd.crosstab(df["Pred"], df["Labels"])

        unique_true = sorted(true_labels.unique())
        subcluster_counts = crosstab.astype(bool).sum(axis=0)
        largest_subclusters = [crosstab[col].idxmax() for col in unique_true]
        max_cells_in_subcluster = crosstab.max(axis=0)
        total_cells = true_labels.value_counts(sort=False)

        avg_cells_per_subcluster = [
            crosstab[col][crosstab[col] != 0].mean() for col in unique_true
        ]
        subcluster_memberships = [
            (crosstab[col][crosstab[col] != 0].index + 1).tolist()
            for col in unique_true
        ]

        summary_df = pd.DataFrame(
            {
                "Number of predicted sub clusters in functional clusters": subcluster_counts,
                "Largest predicted sub cluster in functional clusters": largest_subclusters,
                "Number of cells in largest predicted subcluster": max_cells_in_subcluster,
                "Number of cells in functional clusters": total_cells,
                "Average of cells in all predicted subclusters": avg_cells_per_subcluster,
                "Ratio largest subcluster to functional cluster": (
                    max_cells_in_subcluster / total_cells
                ),
                "Predicted clusters in Labeled clusters": subcluster_memberships,
            },
            index=pd.Index(unique_true, name="Functional clusters"),
        )

        return crosstab, summary_df
    # Format sample metadata and filter out single-cell clusters.
    def prepare_sample_summary(self,summary_df,labels_series,sample_id,patient_id):
        
        cluster_ids, counts = np.unique(labels_series, return_counts=True)
        multi_cell_mask = counts > 1

        filtered_pred_clusters = summary_df[
            "Number of predicted sub clusters in functional clusters"
        ].values[multi_cell_mask]
        filtered_cells_in_fun = summary_df[
            "Number of cells in functional clusters"
        ].values[multi_cell_mask]

        is_apa = 0 if "SN" in sample_id else 1
        sample_tag = f"{sample_id.replace(' ', '_')}_{patient_id}"

        return pd.DataFrame(
            {
                "number_pred_sub_clust_in_fun": filtered_pred_clusters,
                "num_cells_in_fun": filtered_cells_in_fun,
                "sn_apa": is_apa,
                "patient_id": patient_id,
                "sample": sample_tag,
            }
        )
    # Calculate pairwise correlation distributions within and between clusters.
    def calculate_cluster_correlations(self, data, labels_df):

        corr_matrix = data.corr().to_numpy()
        cell_names = data.columns.tolist()

        within_correlations = []
        between_correlations = []

        unique_clusters = labels_df["Clusters"].unique()
        cluster_indices_map = {
            cluster: labels_df[labels_df["Clusters"] == cluster]["Cells"].values
            for cluster in unique_clusters
        }

        for cluster, cells in cluster_indices_map.items():
            idxs = [cell_names.index(c) for c in cells if c in cell_names]
            if len(idxs) > 1:
                sub_corr = corr_matrix[np.ix_(idxs, idxs)]
                triu_indices = np.triu_indices(len(idxs), k=1)
                within_correlations.extend(sub_corr[triu_indices])

        all_triu_indices = np.triu_indices(len(cell_names), k=1)
        all_cells_correlations = corr_matrix[all_triu_indices].tolist()

        return {
            "within_clust_correlations": within_correlations,
            "all_cells_corr": all_cells_correlations,
            "between_cluster_correlation": between_correlations,
        }

    # Finding how many cells are in the predicted sub-clusters.
    def num_cell_in_pred_df(self, crosstab, labels_series, patient_id, sample):
        
        # Exclude single-cell clusters
        cluster_ids, counts = np.unique(labels_series, return_counts=True)
        multi_cell_clusters = cluster_ids[counts > 1]
        
        # Filter crosstab to exclude single-cell functional clusters and zero counts
        filtered_crosstab = crosstab[multi_cell_clusters]
        stacked = filtered_crosstab.unstack().reset_index()
        stacked.columns = ["fun_clust", "pred_clust", "num_cells_pred_(sub)_in_fun"]
        stacked = stacked[stacked["num_cells_pred_(sub)_in_fun"] > 0].copy()

        count_map = dict(zip(cluster_ids, counts))
        stacked["num_cells_in_fun"] = stacked["fun_clust"].map(count_map)
        stacked["sn_apa"] = 0 if "SN" in sample else 1
        stacked["patient_id"] = patient_id
        
        sam_pat = sample.split()
        sample_tag = f"{sam_pat[0]}_{sam_pat[1]}_{patient_id}" if len(sam_pat) > 1 else f"{sample}_{patient_id}"
        stacked["sample"] = sample_tag
        return stacked[[
            "pred_clust", 
            "fun_clust", 
            "num_cells_in_fun", 
            "num_cells_pred_(sub)_in_fun", 
            "sn_apa", 
            "patient_id", 
            "sample"
        ]]

########## Helper module for generation of cluster data visualizations.
class Visualizer:
    # Plot dendrogram
    @staticmethod
    def plot_dendrogram(dataset_name,dissimilarity_matrix,labels,threshold,linkage_method ="complete"): 
        plt.figure(figsize=(12, 5))
        condensed_dist = squareform(dissimilarity_matrix.to_numpy(), checks=False)
        linkage_matrix = linkage(condensed_dist, linkage_method)
        dendrogram(
            linkage_matrix,
            labels=labels,
            orientation="top",
            color_threshold=threshold,
        )
        plt.axhline(y=threshold, color="grey", linewidth=1, linestyle="dashed")
        plt.title(f"Dendrogram for {dataset_name}")
        plt.xlabel("Cells")
        plt.ylabel("Dissimilarity (1-correlation)")
        plt.tight_layout()
        plt.show()

    # Plot a correlation heatmap.
    @staticmethod
    def plot_clustermap(dataset_name, data):
        plt.figure(figsize=(15, 12))
        sns.heatmap(
            data.corr().round(2),
            cmap="RdBu",
            vmin=-1,
            vmax=1,
            xticklabels=True,
            yticklabels=True,
        )
        plt.title(f"Heatmap for {dataset_name}", fontsize=16)
        plt.tight_layout()
        plt.show()