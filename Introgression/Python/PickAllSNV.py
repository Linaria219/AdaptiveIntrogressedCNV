#!/usr/bin/env python3
"""
Function: Processes VCF.gz files of 22 autosomes, extracts the coordinates of non-INDEL SNV loci, and formats the chromosome numbers.
Input: Files chr1.vcf.gz to chr22.vcf.gz in the specified path.
Output: msea.snv.all file containing formatted chromosome numbers and SNV coordinates (tab-separated).
"""

import gzip
import os

def main():
    input_path = "/public/group_data_2023/xiaoyn/Data_MSEA/msea_total/"
    output_file = "msea.snv.all"
    
    #os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, 'w') as out_f:
        # Process chromosomes 1 to 22
        for chrom_num in range(1, 23):
            vcf_file = os.path.join(input_path, f"chr{chrom_num}.vcf.gz")
            
            if not os.path.exists(vcf_file):
                print(f"Warning: File {vcf_file} does not exist, skipping")
                continue
                
            try:
                with gzip.open(vcf_file, 'rt', encoding='utf-8') as in_f:
                    for line in in_f:
                        if line.startswith('#'):
                            continue
                            
                        parts = line.strip().split('\t')
                        # Ensure the line has enough columns
                        if len(parts) < 8:
                            continue
                            
                        # Check if the INFO column (8th column) does not contain "VT=INDEL"
                        info_field = parts[7]
                        if "VT=INDEL" not in info_field:
                            # Extract chromosome number and position
                            chrom = parts[0]
                            pos = parts[1]
                            # Format chromosome number (ensure it is in chr<number> format)
                            formatted_chrom = f"chr{chrom}" if not chrom.startswith('chr') else chrom
                            # Write to output file (tab-separated)
                            out_f.write(f"{formatted_chrom}\t{pos}\n")
                            
            except Exception as e:
                print(f"Error occurred while processing file {vcf_file}: {str(e)}")
                continue
                
    print(f"Processing complete! Results saved to {output_file}")

if __name__ == "__main__":
    main()
