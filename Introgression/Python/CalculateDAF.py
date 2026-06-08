#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DAF (derived allele frequency) calculation script
Function: Calculates DAF values for target, AFR, and ARC populations based on VCF files.
Inputs: VCF file, AA reference file, ARC file, sample list
Output: File containing DAF values for the three populations.
"""

import gzip
import os
import sys
import logging
import bisect
from collections import defaultdict
from typing import Dict, List, Tuple, Set, Optional

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def read_sample_list(file_path: str) -> Set[str]:
    """Read sample list file"""
    samples = set()
    try:
        with open(file_path, 'r') as f:
            for line in f:
                sample = line.strip()
                if sample:
                    samples.add(sample)
        logger.info(f"Read {len(samples)} samples from file {file_path}")
        return samples
    except Exception as e:
        logger.error(f"Error reading sample list file {file_path}: {str(e)}")
        raise

def read_aa_reference(aa_bed_path: str) -> Dict[int, str]:
    """Read AA(ancestral allele) reference file"""
    aa_dict = {}
    if not os.path.exists(aa_bed_path):
        logger.warning(f"AA reference file does not exist: {aa_bed_path}")
        return aa_dict
    
    try:
        with open(aa_bed_path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 4:
                    chrom = fields[0].replace('chr', '')
                    pos = int(fields[2])
                    aa_base = fields[3].upper()
                    aa_dict[pos] = aa_base
        
        logger.info(f"Loaded {len(aa_dict)} AA status records from file {aa_bed_path}")
        return aa_dict
    except Exception as e:
        logger.error(f"Error reading AA reference file {aa_bed_path}: {str(e)}")
        raise

def calculate_daf(genotypes: List[str], ref: str, alt: str, aa_base: str) -> Optional[float]:
    """Calculate DAF value"""
    # Normalize base case
    ref = ref.upper()
    alt = alt.upper()
    aa_base = aa_base.upper()
    
    # Identify the non-AA allele
    if ref == aa_base:
        non_aa_allele = alt
    elif alt == aa_base:
        non_aa_allele = ref
    else:
        # Check if ALT contains multiple alleles
        alt_alleles = alt.split(',')
        if aa_base in alt_alleles:
            non_aa_allele = ref
        else:
            return None  # AA is not in REF or ALT
    
    # Count the non-AA alleles
    non_aa_count = 0
    total_alleles = 0
    
    for gt in genotypes:
        # Normalize genotype format
        if '|' in gt:
            alleles = gt.split('|')
        elif '/' in gt:
            alleles = gt.split('/')
        else:
            continue
        
        # Count alleles
        for allele in alleles:
            if allele == '0':  # REF allele
                if non_aa_allele == ref:
                    non_aa_count += 1
            elif allele == '1':  # First ALT allele
                if non_aa_allele == alt.split(',')[0]:
                    non_aa_count += 1
            # Ignore other alleles (only process biallelic sites)
        
        total_alleles += 2
    
    if total_alleles == 0:
        return None
    
    daf = non_aa_count / total_alleles
    return daf

def process_target_vcf_with_aa(chrom: str, vcf_path: str, aa_dict: Dict[int, str], 
                            target_samples: Set[str], output_dir: str) -> str:
    """Process target VCF file and calculate DAF, generate intermediate file 1"""
    # Prepare output file
    output_file = os.path.join(output_dir, f"temp.{chrom}.AA.target.deni.daf")
    
    try:
        # Determine file opening method
        open_func = gzip.open if vcf_path.endswith('.gz') else open
        mode = 'rt' if vcf_path.endswith('.gz') else 'r'
        
        with open_func(vcf_path, mode) as vcf_file, open(output_file, 'w') as out_f:
            # Read header line, find sample column indices
            sample_indices = []
            sample_names = []
            
            for line in vcf_file:
                if line.startswith('#CHROM'):
                    header_fields = line.strip().split('\t')
                    sample_names = header_fields[9:]
                    
                    # Find the column indices for target samples
                    for i, sample_name in enumerate(sample_names):
                        if sample_name in target_samples:
                            sample_indices.append(i)
                    
                    logger.info(f"Find {len(sample_indices)} target samples")
                    break
            
            # Process data rows
            processed_count = 0
            valid_count = 0
            
            for line in vcf_file:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 10:
                    continue
                
                # Get basic information
                chrom_num = fields[0].replace('chr', '')
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                
                # Check if the position has AA information
                if pos not in aa_dict:
                    continue
                
                aa_base = aa_dict[pos]
                
                # Extract target sample genotypes
                target_genotypes = []
                for i in sample_indices:
                    if i < len(fields) - 9:
                        gt_field = fields[9 + i].split(':')[0]
                        target_genotypes.append(gt_field)
                
                # Calculate DAF
                daf = calculate_daf(target_genotypes, ref, alt, aa_base)
                if daf is None:
                    continue
                
                # Write to intermediate file
                out_f.write(f"{pos}\t{aa_base}\t{daf:.6f}\n")
                valid_count += 1
                processed_count += 1
                
                if processed_count % 100000 == 0:
                    logger.info(f"Processed {processed_count} target data rows, valid sites: {valid_count}")
            
            logger.info(f"Completed target processing: Total sites processed: {processed_count}, Valid sites: {valid_count}")
            return output_file
    
    except Exception as e:
        logger.error(f"Error occurred while processing target VCF file: {str(e)}")
        raise

def process_afr_vcf_with_intermediate(chrom: str, vcf_path: str, intermediate_file: str, 
                                     afr_samples: Set[str], output_dir: str) -> str:
    """Process AFR VCF file and match with intermediate file 1, generate intermediate file 2"""
    # Prepare output file
    output_file = os.path.join(output_dir, f"temp.{chrom}.AA.target.afr.deni.daf")
    
    # Read intermediate file 1 data
    intermediate_data = {}
    with open(intermediate_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            pos = int(fields[0])
            aa_base = fields[1]
            target_daf = float(fields[2])
            intermediate_data[pos] = (aa_base, target_daf)
    
    logger.info(f"Loaded {len(intermediate_data)} valid sites from intermediate file 1")
    
    try:
        # Determine file opening method
        open_func = gzip.open if vcf_path.endswith('.gz') else open
        mode = 'rt' if vcf_path.endswith('.gz') else 'r'
        
        with open_func(vcf_path, mode) as vcf_file, open(output_file, 'w') as out_f:
            # Read header line and find sample column indices
            sample_indices = []
            sample_names = []
            
            for line in vcf_file:
                if line.startswith('#CHROM'):
                    header_fields = line.strip().split('\t')
                    sample_names = header_fields[9:]
                    
                    # Find column indices for AFR samples
                    for i, sample_name in enumerate(sample_names):
                        if sample_name in afr_samples:
                            sample_indices.append(i)
                    
                    logger.info(f"Found {len(sample_indices)} AFR samples")
                    break
            
            # Process data rows
            processed_count = 0
            valid_count = 0
            intermediate_positions = sorted(intermediate_data.keys())
            pos_index = 0
            
            for line in vcf_file:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 10:
                    continue
                
                # Get basic information
                chrom_num = fields[0].replace('chr', '')
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                
                # Skip non-SNP sites
                if 'VT=SNP' not in fields[7]:
                    continue
                
                # Skip positions not in the intermediate file
                while pos_index < len(intermediate_positions) and intermediate_positions[pos_index] < pos:
                    pos_index += 1
                
                if pos_index >= len(intermediate_positions) or intermediate_positions[pos_index] != pos:
                    continue
                
                # Get AA information
                aa_base, target_daf = intermediate_data[pos]
                
                # Extract AFR sample genotypes
                afr_genotypes = []
                for i in sample_indices:
                    if i < len(fields) - 9:
                        gt_field = fields[9 + i].split(':')[0]
                        afr_genotypes.append(gt_field)
                
                # Calculate DAF
                daf = calculate_daf(afr_genotypes, ref, alt, aa_base)
                if daf is None:
                    continue
                
                # Write to intermediate file 2
                out_f.write(f"{pos}\t{aa_base}\t{target_daf}\t{daf:.6f}\n")
                valid_count += 1
                processed_count += 1
                
                if processed_count % 100000 == 0:
                    logger.info(f"Processed {processed_count} lines of AFR data, valid sites: {valid_count}")
            
            logger.info(f"Completed AFR processing: Total sites processed: {processed_count}, Valid sites: {valid_count}")
            return output_file
    
    except Exception as e:
        logger.error(f"Error occurred while processing AFR VCF file: {str(e)}")
        raise

def process_arc_with_intermediate(chrom: str, arc_path: str, intermediate_file: str, output_dir: str) -> str:
    """Process ARC file and match with intermediate file 2 to generate final file"""
    # Prepare output file
    output_file = os.path.join(output_dir, f"chr{chrom}.target.afr.deni.daf")
    
    # Read data from intermediate file 2
    intermediate_data = {}
    with open(intermediate_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue
            pos = int(fields[0])
            aa_base = fields[1]
            target_daf = float(fields[2])
            afr_daf = float(fields[3])
            intermediate_data[pos] = (aa_base, target_daf, afr_daf)
    
    logger.info(f"Loaded {len(intermediate_data)} valid sites from intermediate file 2")
    
    try:
        # Determine file opening method
        open_func = gzip.open if arc_path.endswith('.gz') else open
        mode = 'rt' if arc_path.endswith('.gz') else 'r'
        
        with open_func(arc_path, mode) as arc_file, open(output_file, 'w') as out_f:
            # Write header
            out_f.write("#CHROM\tPOS\ttarget_DAF\tAFR_DAF\tARC_DAF\n")
            
            # Skip header lines
            for line in arc_file:
                if line.startswith('#CHROM'):
                    break
            
            # Process data lines
            processed_count = 0
            valid_count = 0
            intermediate_positions = sorted(intermediate_data.keys())
            pos_index = 0
            
            for line in arc_file:
                fields = line.strip().split('\t')
                if len(fields) < 10:
                    continue
                
                # Get basic information
                chrom_num = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                
                # Skip positions not in the intermediate file
                while pos_index < len(intermediate_positions) and intermediate_positions[pos_index] < pos:
                    pos_index += 1
                
                if pos_index >= len(intermediate_positions) or intermediate_positions[pos_index] != pos:
                    continue
                
                # Get AA information
                aa_base, target_daf, afr_daf = intermediate_data[pos]
                
                # Extract ARC genotype
                gt_field = fields[9].split(':')[0]
                
                # Calculate ARC DAF
                arc_daf = calculate_arc_daf(gt_field, ref, alt, aa_base)
                if arc_daf is None:
                    continue
                
                # Write to final file
                out_f.write(f"chr{chrom}\t{pos}\t{target_daf:.6f}\t{afr_daf:.6f}\t{arc_daf:.6f}\n")
                valid_count += 1
                processed_count += 1
                
                if processed_count % 100000 == 0:
                    logger.info(f"Processed {processed_count} lines of ARC data, Valid sites: {valid_count}")
            
            logger.info(f"Completed ARC processing: Total sites processed: {processed_count}, Valid sites: {valid_count}")
            return output_file
    
    except Exception as e:
        logger.error(f"Error occurred while processing ARC file: {str(e)}")
        raise

def calculate_arc_daf(gt: str, ref: str, alt: str, aa_base: str) -> Optional[float]:
    """Calculate ARC DAF value (single individual)"""
    # Unified base case
    ref = ref.upper()
    alt = alt.upper()
    aa_base = aa_base.upper()
    
    # Determine non-AA allele
    if ref == aa_base:
        non_aa_allele = alt
    elif alt == aa_base:
        non_aa_allele = ref
    else:
        # Check if ALT contains multiple alleles
        alt_alleles = alt.split(',')
        if aa_base in alt_alleles:
            non_aa_allele = ref
        else:
            return None  # AA not in REF or ALT
    
    # Standardize genotype format
    if '|' in gt:
        alleles = gt.split('|')
    elif '/' in gt:
        alleles = gt.split('/')
    else:
        return None
    
    # Calculate non-AA allele count
    non_aa_count = 0
    
    for allele in alleles:
        if allele == '0':  # REF allele
            if non_aa_allele == ref:
                non_aa_count += 1
        elif allele == '1':  # First ALT allele
            if non_aa_allele == alt.split(',')[0]:
                non_aa_count += 1
        # Ignore other alleles (only process biallelic)
    
    # Calculate DAF (single individual)
    daf = non_aa_count / 2.0
    return daf

def main():
    """Main function"""
    # Path configuration
    base_dir = "Path/to/your/data"
    output_dir = "Path/to/your/output"
    temp_dir = os.path.join(output_dir, "temp")
    
    # Ensure directory exists
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    
    # Sample list file
    target_sample_file = os.path.join(base_dir, "sample.target.list")
    afr_sample_file = os.path.join(base_dir, "sample.afr.list")
    
    # Read sample lists
    target_samples = read_sample_list(target_sample_file)
    afr_samples = read_sample_list(afr_sample_file)
    
    # Process all autosomes (1-22)
    chromosomes = [str(i) for i in range(1, 9)]
    
    for chrom in chromosomes:
        logger.info(f"Starting processing chromosome {chrom}")
        
        # File paths
        aa_bed_path = f"path/to/ancestor/{chrom}ancestor.bed"
        target_vcf_path = f"path/to/targetVCF/chr{chrom}.phased.vcf.gz"
        afr_vcf_path = f"path/to/AFRVCF/chrom_{chrom}.vcf.gz"
        arc_path = f"path/to/ARC/chrom_{chrom}.vcf.gz"
        
        # Step 1: Read AA reference file
        aa_dict = read_aa_reference(aa_bed_path)
        if not aa_dict:
            logger.warning(f"Chromosome {chrom} has no valid AA data, skipping processing")
            continue
        
        # Step 2: Process target VCF and generate intermediate file 1
        logger.info("Step 1: Processing target VCF")
        intermediate1 = process_target_vcf_with_aa(chrom, target_vcf_path, aa_dict, target_samples, temp_dir)
        
        # Step 3: Process AFR VCF and generate intermediate file 2
        logger.info("Step 2: Processing AFR VCF")
        intermediate2 = process_afr_vcf_with_intermediate(chrom, afr_vcf_path, intermediate1, afr_samples, temp_dir)
        
        # Step 4: Process ARC file and generate final file
        logger.info("Step 3: Processing ARC file")
        final_output = process_arc_with_intermediate(chrom, arc_path, intermediate2, output_dir)
        
        logger.info(f"Completed chromosome {chrom}: Final results saved in {final_output}")
        
        # Clean up intermediate files (optional)
        # os.remove(intermediate1)
        # os.remove(intermediate2)
    
    logger.info("All chromosomes processed!")

if __name__ == "__main__":
    main()
