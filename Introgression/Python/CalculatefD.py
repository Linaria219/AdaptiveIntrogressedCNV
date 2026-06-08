#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fD statistic calculation script. 
Function: Calculates the fD statistic based on the DAF value of the target, AFR(africa), and ARC(archaic) populations.
"""

import os
import gzip
import logging
import numpy as np
from typing import List, Tuple, Optional

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def read_daf_file(file_path: str) -> List[Tuple[str, int, float, float, float]]:
    """
    Read and parse the DAF(derived allele frequency) file to extract valid SNV data
    
    Args:
        file_path: Path to the DAF file
        
    Returns:
        A list of valid SNV data, each element is a tuple of (chrom, pos, target_daf, afr_daf, arc_daf)
    """
    data = []
    try:
        if not os.path.exists(file_path):
            logger.error(f"DAF file does not exist: {file_path}")
            return data
        
        # Determine the file opening method
        if file_path.endswith('.gz'):
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'
        
        with open_func(file_path, mode) as f:
            for line in f:
                if not line.startswith('#'):
                    break
           
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 5:
                    continue
                
                try:
                    chrom = fields[0]
                    pos = int(fields[1])
                    target_daf = fields[2]
                    afr_daf = fields[3]
                    arc_daf = fields[4]
                    
                    # Only keep SNVs where all DAF values are valid
                    if (target_daf != '.' and afr_daf != '.' and arc_daf != '.'):
                        target_val = float(target_daf)
                        afr_val = float(afr_daf)
                        arc_val = float(arc_daf)
                        
                        # Validate DAF values are within a reasonable range [0, 1]
                        if (0 <= target_val <= 1 and 
                            0 <= afr_val <= 1 and 
                            0 <= arc_val <= 1):
                            data.append((chrom, pos, target_val, afr_val, arc_val))
                
                except ValueError as e:
                    logger.warning(f"Ignoring line {line_num}: Unable to parse numeric value - {str(e)}")
                    continue
        
        logger.info(f"{len(data)} valid SNVs were read from file {file_path}.")
        return data
        
    except Exception as e:
        logger.error(f"Error reading DAF file {file_path}: {str(e)}")
        raise

def calculate_fd_window(window_data: List[Tuple[float, float, float]]) -> Optional[float]:
    """
    Calculate the fD statistic for a given window of SNV data
    
    Args:
        window_data: A list of DAF values for each SNV in the window, each element is a tuple of (target_daf, afr_daf, arc_daf)
        
    Returns:
        fD value or None (if unable to calculate)
    """
    if not window_data:
        return None
    
    numerator_sum = 0.0
    denominator_sum = 0.0
    
    for target_daf, afr_daf, arc_daf in window_data:
        # Calculate Pd (maximum value between target and arc)
        pd = max(target_daf, arc_daf)
        
        # Calculate numerator: (1-AFR)*target*ARC - AFR*(1-target)*ARC
        numerator = ((1 - afr_daf) * target_daf * arc_daf - afr_daf * (1 - target_daf) * arc_daf)
        
        # Calculate denominator: (1-AFR)*Pd*Pd - AFR*(1-Pd)*Pd
        denominator = ((1 - afr_daf) * pd * pd - afr_daf * (1 - pd) * pd)
        
        numerator_sum += numerator
        denominator_sum += denominator
    
    # Avoid division by zero
    if abs(denominator_sum) < 1e-10:
        return None
    
    fd = numerator_sum / denominator_sum
    return fd

def process_chromosome(chrom: str, input_dir: str, output_dir: str, window_size: int = 100, step_size: int = 50) -> None:
    """
    Process a single chromosome's DAF file and calculate fD values
    
    Args:
        chrom: Chromosome number
        input_dir: Input file directory
        output_dir: Output file directory
        window_size: Window size
        step_size: Sliding step size
    """
    # Construct file path
    input_filename = f"chr{chrom}.target.afr.arc.daf"
    input_path = os.path.join(input_dir, input_filename)
    
    output_filename = f"chr{chrom}.target.arc.fd"
    output_path = os.path.join(output_dir, output_filename)
    
    logger.info(f"Starting to process chromosome {chrom}")
    
    # Read DAF file
    snv_data = read_daf_file(input_path)
    if not snv_data:
        logger.warning(f"Chromosome {chrom} has no valid data, skipping processing")
        return
    
    # Prepare output file
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Use sliding window to calculate fD
    results = []
    total_snvs = len(snv_data)
    
    for start_idx in range(0, total_snvs - window_size + 1, step_size):
        end_idx = start_idx + window_size
        
        # Extract window data
        window_snvs = snv_data[start_idx:end_idx]
        window_data = [(target, afr, arc) for _, _, target, afr, arc in window_snvs]
        
        # Calculate fD value
        fd_value = calculate_fd_window(window_data)
        
        if fd_value is not None:
            # Get the coordinate range of the window
            min_pos = window_snvs[0][1] - 1  # Minimum coordinate - 1
            max_pos = window_snvs[-1][1]     # Maximum coordinate
            
            results.append(f"{chrom}\t{min_pos}\t{max_pos}\t{fd_value:.6f}\n")
    
    # Write results
    if results:
        with open(output_path, 'w') as out_file:
            out_file.write("#CHROM\tSTART\tEND\tfD\n")
            out_file.writelines(results)
        
        logger.info(f"Completed chromosome {chrom}: Calculated fD values for {len(results)} windows")
    else:
        logger.warning(f"Chromosome {chrom}: No valid fD values calculated for any window")

def main():
    """Main function"""
    # Path settings
    input_base_dir = "path/to/your/data"
    output_base_dir = "path/to/your/output"
    
    # Process all chromosomes (1-22 and X, Y, etc., adjust as needed)
    chromosomes = [str(i) for i in range(1, 23)]  # Autosomes
    # If you need to process sex chromosomes, you can add them:
    # chromosomes += ['X', 'Y']
    
    for chrom in chromosomes:
        process_chromosome(chrom, input_base_dir, output_base_dir)
    
    logger.info("All chromosomes processed!")

if __name__ == "__main__":
    main()
