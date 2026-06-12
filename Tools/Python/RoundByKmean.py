#!/usr/bin/env python3
"""
CNV Data Adaptive K-Means Rounding Script - A New Logic Based on Central Value Difference Merging
"""

import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from multiprocessing import Pool, cpu_count
import os
import sys
from typing import List, Tuple

def is_integer_array(arr, tol=1e-6):
    """Check if the array contains only integers (vectorized operation)"""
    return np.allclose(arr, np.round(arr), atol=tol)

def adaptive_center_merging_rounding(row_data, min_center_diff=0.5):
    """
    A new rounding logic based on merging cluster centers by their differences
    Input: row_data (1D array, e.g., CN values from 1889 samples)
    Output: rounded_cn (1D array after rounding)
    """
    values = np.array(row_data, dtype=float)
    
    # Check if all values are already integers
    if is_integer_array(values):
        return values.tolist()
    
    # Remove NaN values
    non_nan_mask = ~np.isnan(values)
    valid_values = values[non_nan_mask]
    
    if len(valid_values) <= 1:
        return row_data  # Data insufficient, return original values
    
    # 1. Determine CN range and initial K value for the new logic
    min_val = np.min(valid_values)
    max_val = np.max(valid_values)
    
    # Floor the minimum value and ceil the maximum value
    cn_min = int(np.floor(min_val))
    cn_max = int(np.ceil(max_val))
    initial_k = max(1, cn_max - cn_min + 1)  # Ensure at least 1 cluster
    
    # Limit the maximum K value to avoid over-segmentation
    max_allowed_k = min(10, len(valid_values) - 1)
    initial_k = min(initial_k, max_allowed_k)
    
    if initial_k == 1:
        # Single cluster, round directly
        rounded_valid = np.round(valid_values)
    else:
        # 2. Use the initial K value for K-Means clustering
        kmeans = KMeans(n_clusters=initial_k, random_state=42, n_init=10)
        X = valid_values.reshape(-1, 1)
        
        try:
            clusters = kmeans.fit_predict(X)
            centers = kmeans.cluster_centers_.flatten()
            
            # 3. Central value merging logic
            optimal_centers = merge_centers_by_difference(centers, min_diff=min_center_diff)
            final_k = len(optimal_centers)
            
            if final_k == 1:
                # After merging, only one center remains, so we directly round down.
                rounded_valid = np.round(valid_values)
            else:
                # 4. Use the merged centers to reassign samples
                rounded_valid = assign_samples_to_centers(valid_values, optimal_centers)
                
        except Exception as e:
            # Fallback to simple rounding when clustering fails
            rounded_valid = np.round(valid_values)
    
    # Place the results back in their original positions (keeping NaN positions)
    result = np.full_like(values, np.nan)
    result[non_nan_mask] = rounded_valid
    
    return result.tolist()

def merge_centers_by_difference(centers: np.ndarray, min_diff: float = 0.5) -> List[float]:
    """
    Merge adjacent cluster centers whose difference is less than the threshold
    Return the list of merged center values
    """
    if len(centers) <= 1:
        return centers.tolist()
    
    # Sort the center values
    sorted_centers = np.sort(centers)
    merged_centers = []
    
    i = 0
    while i < len(sorted_centers):
        if i == len(sorted_centers) - 1:
            # Last center, add directly
            merged_centers.append(sorted_centers[i])
            break
        
        # Check the difference between the current center and the next center
        current_center = sorted_centers[i]
        next_center = sorted_centers[i + 1]
        difference = next_center - current_center
        
        if difference < min_diff:
            # Difference is too small, merge these two centers
            merged_center = (current_center + next_center) / 2.0
            merged_centers.append(merged_center)
            i += 2  # Skip the next center, as it has been merged
        else:
            # Difference is large enough, keep the current center
            merged_centers.append(current_center)
            i += 1
    
    # Recursively check, until no more centers can be merged
    if len(merged_centers) < len(centers):
        return merge_centers_by_difference(np.array(merged_centers), min_diff)
    else:
        return merged_centers

def assign_samples_to_centers(samples: np.ndarray, centers: List[float]) -> np.ndarray:
    """
    Assign samples to the nearest center and then round the center values
    """
    centers_array = np.array(centers)
    rounded_centers = np.round(centers_array)  # Round the center values
    
    # Calculate the distance from each sample to each rounded center
    distances = np.abs(samples[:, np.newaxis] - rounded_centers)
    
    # Find the index of the nearest center
    nearest_center_indices = np.argmin(distances, axis=1)
    
    # Assign the rounded center values
    return rounded_centers[nearest_center_indices]

def process_single_row_wrapper(args):
    """
    Wrapper function for multiprocessing
    Input: (row_index, row_data, header_info)
    Output: (row_index, processed_result, error_msg)
    """
    row_index, row_data, header_info = args
    try:
        # Extract coordinate information (first 3 columns)
        coord_info = row_data[:3].tolist()
        # Extract CN values (4th column to the end)
        cn_values = row_data[3:].tolist()
        
        # Use the new rounding logic to process CN values
        rounded_cn = adaptive_center_merging_rounding(cn_values, min_center_diff=0.5)
        
        # Combine the results
        result_row = coord_info + rounded_cn
        return (row_index, result_row, None)
    except Exception as e:
        return (row_index, None, str(e))

def main_optimized(input_file, output_file, n_jobs=10, chunk_size=100):
    """
    Main optimized function: new biological logic + multiprocessing + batch I/O
    """
    print(f"Starting to read file: {input_file}")
    print(f"Using center difference merging algorithm (minimum center difference=0.5)")
    
    try:
        # Read the file
        df = pd.read_csv(input_file, sep='\t', header=0)
    except Exception as e:
        print(f"Failed to read file {input_file}: {e}")
        return False
    
    if df.shape[1] < 4:
        print("Error: File has fewer than 4 columns, which does not meet the requirements")
        return False
    
    print(f"File contains {df.shape[0]} rows and {df.shape[1]} columns")
    print(f"Using {n_jobs} processes for parallel computation, batch size is {chunk_size}")
    
    # Prepare multiprocessing parameters
    rows_data = [(i, row, df.columns.tolist()) for i, (idx, row) in enumerate(df.iterrows())]
    
    # Create output file (write header first)
    df.iloc[:0].to_csv(output_file, sep='\t', index=False)
    
    # Process in batches to reduce memory usage
    total_batches = (len(rows_data) + chunk_size - 1) // chunk_size
    processed_count = 0
    
    for batch_num in range(total_batches):
        start_idx = batch_num * chunk_size
        end_idx = min((batch_num + 1) * chunk_size, len(rows_data))
        batch_data = rows_data[start_idx:end_idx]
        
        print(f"Processing batch {batch_num + 1}/{total_batches} (rows {start_idx + 1}-{end_idx})")
        
        # Use multiprocessing to process the current batch
        with Pool(processes=n_jobs) as pool:
            results = pool.map(process_single_row_wrapper, batch_data)
        
        # Process the results
        batch_results = []
        errors = []
        
        for row_index, result_row, error_msg in results:
            if error_msg is None and result_row is not None:
                batch_results.append((row_index, result_row))
            else:
                errors.append((row_index, error_msg))
                # Keep the original line when an error occurs
                original_row = df.iloc[row_index].tolist()
                batch_results.append((row_index, original_row))
        
        # Sort by original order and extract results
        batch_results.sort(key=lambda x: x[0])
        ordered_results = [result for _, result in batch_results]
        
        # Batch write the current batch results
        result_batch_df = pd.DataFrame(ordered_results, columns=df.columns)
        result_batch_df.to_csv(output_file, mode='a', sep='\t', 
                             header=False, index=False)
        
        # Report errors
        if errors:
            print(f"Batch {batch_num + 1} found {len(errors)} errors")
        
        processed_count += len(batch_data)
    
    print(f"Processing complete! Total rows processed: {processed_count}")
    print(f"Results saved to: {output_file}")
    return True

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python optimized_cnv_rounding.py <input_file> [n_jobs] [chunk_size]")
        print("Parameters:")
        print("  <input_file> : Input CNV data file (TAB separated)")
        print("  [n_jobs]     : Number of parallel processes (default: 10)")
        print("  [chunk_size] : Batch processing size (default: 100)")
        sys.exit(1)
    
    input_file = sys.argv[1]
    n_jobs = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    chunk_size = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    
    # Generate output file name
    base_name = os.path.splitext(input_file)[0]
    output_file = f"{base_name}_kmean_optimized.bed"
    
    # Run optimized processing
    success = main_optimized(input_file, output_file, n_jobs, chunk_size)
    
    if success:
        print("Processing successful!")
    else:
        print("Processing failed. Please check the error messages.")
