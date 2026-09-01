from pathlib import Path
from typing import Dict, Any
import pandas as pd
import numpy as np
from scipy.signal import savgol_filter
import pickle
from patient_clustering_analyzer import PatientAnalyzer, Visualizer

DATA_DIR = Path("./data")
OUTPUT_DIR = Path("./output")
# Cutt off distance for hierarchical clustering  
DISTANCE_THRESHOLD = 0.15 
# Filtering method
FILTER_METHOD = "savgol_filter"
# End of valid time series in (s)
CUTOFF_CONFIG = {
    "Patient_38": {"SN 3": 319, "SN 5": 274, "APA 4": 355, "APA 6": 314},
    "Patient_41": {
        "SN 1": 669,
        "SN 3": 499,
        "SN 5": 899,
        "APA 2": 709,
        "APA 4": 599,
    },
    "Patient_42": {
        "SN 1": 553,
        "SN 4": 464,
        "SN 6": 504,
        "SN 7": 869,
        "APA 2": 679,
        "APA 3": 419,
        "APA 5": 429,
    },
    "Patient_43": {"SN 2": 1159, "SN 4": 984},
    "Patient_45": {"SN 2": 899, "APA 4": 739},
    "Patient_46": {"SN 4": 479, "SN 5": 709, "APA 2": 879},
    "Patient_48": {"SN 1": 839, "SN 3": 539, "APA 2": 724, "APA 4": 626},
    "Patient_49": {"APA 1": 689},
    "Patient_53": {"SN 1": 619, "APA 2": 507},
}
# Differnet filter params
FILTER_CONFIG = {
    "mean_filter":{"window":15,"center":True},
    "savgol_filter":{"window":13,"polyorder":3}    
}

########## Normalize and filter raw time-series trace signals

# df: Raw input signal dataframe.
# cutoff: End of valid time series in (s). (cuttoff_config)
# filter_type: Type of digital filter applied (savgol_filter or mean_filter).

def preprocess_signal_data(df, cutoff, filter_type):

    sliced_df = df.iloc[:, 1:].loc[df.index <= cutoff]
    normalized = (sliced_df - sliced_df.mean()) / sliced_df.std()

    if filter_type == "mean_filter":
        return normalized.rolling(window=FILTER_CONFIG["mean_filter"]["window"], center=FILTER_CONFIG["mean_filter"]["center"]).mean()
    elif filter_type == "savgol_filter":
        filtered_array = savgol_filter(
            normalized, window_length=FILTER_CONFIG["mean_filter"]["window"], polyorder=FILTER_CONFIG["mean_filter"]["polyorder"], axis=0)
        return pd.DataFrame(
            filtered_array, index=sliced_df.index, columns=sliced_df.columns)
    return normalized

#################### Function to run pipeline on all patients 
def run_pipeline():
    #Custom class module for analyzing patient single-cell calcium trace clustering data
    analyzer = PatientAnalyzer()

    processed_results = {}

    for patient_id, samples in CUTOFF_CONFIG.items():
        file_path = DATA_DIR / f"{patient_id}.xlsx"
        labels_path = DATA_DIR / f"{patient_id}_labels.xlsx"

        patient_xls = pd.read_excel(file_path, sheet_name=None)
        labels_xls = pd.read_excel(labels_path, sheet_name=None)
        processed_results[patient_id] = {}

        for sample_name, raw_data in patient_xls.items():
            cutoff = samples[sample_name]
            clean_data = preprocess_signal_data(raw_data, cutoff, FILTER_METHOD)
            labels = labels_xls[sample_name]

            # Clustering analysis
            predictions = analyzer.predict_clusters(clean_data, DISTANCE_THRESHOLD)
            ct, summary_df = analyzer.analyze_clusters(
                predictions, labels["Clusters"]
            )
            sample_summary = analyzer.prepare_sample_summary(
                summary_df, labels["Clusters"], sample_name, patient_id
            )
            correlations = analyzer.calculate_cluster_correlations(clean_data, labels)
            number_of_cells_in_predicted_clusters = analyzer.num_cell_in_pred_df(ct, labels, patient_id, sample_name)
            processed_results[patient_id][sample_name] = {
                "predictions": predictions,
                "summary": sample_summary,
                "correlations": correlations,
                "number_of_cells_in_predicted_clusters":number_of_cells_in_predicted_clusters 
            }
            print(f"Processed {patient_id} - {sample_name}")
    processed_results["config"] = {"filter_method":FILTER_METHOD,"filter_config":FILTER_CONFIG[FILTER_METHOD],"Distance":DISTANCE_THRESHOLD}        
    return processed_results


if __name__ == "__main__":
    results = run_pipeline()
    with open(f"clustering_filter_{FILTER_METHOD}_distance_{DISTANCE_THRESHOLD}.pkl","wb") as f:
        pickle.dump(results,f)