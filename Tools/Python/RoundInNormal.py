#!/usr/bin/env python3
from decimal import Decimal, ROUND_HALF_UP
import sys
import os

def traditional_round(x):
    """
    Implement traditional rounding (rounding up when the number reaches five)
    Use the Decimal module to avoid floating-point precision issues.
    """
    try:
        d = Decimal(str(x))
        rounded = d.quantize(Decimal('1'), rounding=ROUND_HALF_UP)
        return int(rounded)
    except:
        return None  # Conversion failed, return None

def is_numeric(value):
    """
    Check if the value is a valid number
    """
    try:
        float(value)
        return True
    except (ValueError, TypeError):
        return False

def validate_row_data(cn_values_str):
    """
    Check if all CN values in a row are valid numbers
    """
    for val in cn_values_str:
        if not is_numeric(val):
            return False
    return True

def determine_cnv_type(cn_values):
    """
    Determine the CNV type based on CN values
    """
    if all(cn == 2 for cn in cn_values):
        return "norm"
    elif all(cn <= 2 for cn in cn_values):
        return "DEL"
    elif all(2 <= cn <= 4 for cn in cn_values):
        return "DUP"
    else:
        return "mCNV"

def main():
    input_file = "path/to/Old.bed"
    output_file = "path/to/Processed.bed"
    
    try:
        with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
            # Write the header for the output file
            f_out.write("CHROM\tSTART\tEND\tTYPE\tNUM_TYPES\t" + 
                       "\t".join([f"CN_{i+1}" for i in range(1000)]) + "\n")
            
            line_count = 0
            processed_count = 0
            skipped_count = 0
            
            for line in f_in:
                line = line.strip()
                if not line:
                    continue
                    
                # Process comment lines (lines starting with #)
                if line.startswith('#'):
                    f_out.write(line + '\n')
                    continue
                    
                parts = line.split('\t')
                if len(parts) < 4:
                    print(f"Warning: Line {line_count+1} has insufficient columns, skipping")
                    skipped_count += 1
                    continue
                
                # Extract the first three columns (chromosome, start, end)
                chrom, start, end = parts[0], parts[1], parts[2]
                
                # Extract CN value columns (starting from the 4th column)
                cn_values_str = parts[3:]
                
                # Quality control: Check if all CN values in the row are valid numbers
                if not validate_row_data(cn_values_str):
                    print(f"Warning: Line {line_count+1} contains non-numeric content, skipping")
                    skipped_count += 1
                    line_count += 1
                    continue
                
                # Perform traditional rounding on each CN value
                rounded_ints = []
                valid_conversion = True
                
                for val in cn_values_str:
                    rounded_val = traditional_round(val)
                    if rounded_val is None:
                        print(f"Error: Line {line_count+1} value '{val}' conversion failed")
                        valid_conversion = False
                        break
                    rounded_ints.append(rounded_val)
                
                if not valid_conversion:
                    print(f"Warning: Line {line_count+1} contains invalid CN values, skipping")
                    skipped_count += 1
                    line_count += 1
                    continue
                
                # Statistics: Count unique CN values
                unique_cn = set(rounded_ints)
                num_types = len(unique_cn)
                
                # Determine CNV type
                cnv_type = determine_cnv_type(rounded_ints)
                
                rounded_str = [str(x) for x in rounded_ints]
                # Only output actual columns, avoiding excessive empty columns
                output_parts = [chrom, start, end, cnv_type, str(num_types)] + rounded_str
                f_out.write('\t'.join(output_parts) + '\n')
                
                processed_count += 1
                line_count += 1
                
                if line_count % 1000 == 0:
                    print(f"Processed {line_count} lines, succeeded {processed_count} lines, skipped {skipped_count} lines")
            
            print(f"Processing complete! Total lines processed: {line_count}")
            print(f"Successful processing: {processed_count} lines")
            print(f"Skipped processing: {skipped_count} lines (containing non-numeric content or conversion errors)")
            print(f"Results saved to: {output_file}")
            
    except FileNotFoundError:
        print(f"Error: Input file does not exist: {input_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing file {input_file}: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
