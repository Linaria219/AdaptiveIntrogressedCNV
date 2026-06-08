import pandas as pd
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
import argparse

def fill_zero_k_values_by_gc_direction(k_values, bins, gc_threshold_low=0.25, gc_threshold_high=0.75):
    k_values = k_values.copy()
    n = len(k_values)
    for i in range(n):
        if np.isnan(k_values[i]) or k_values[i] == 0:
            gc_mid = (bins[i] + bins[i+1]) / 2
            if gc_mid < gc_threshold_low:
                for j in range(i+1, n):
                    if not np.isnan(k_values[j]) and k_values[j] != 0:
                        k_values[i] = k_values[j]
                        break
            elif gc_mid > gc_threshold_high:
                for j in range(i-1, -1, -1):
                    if not np.isnan(k_values[j]) and k_values[j] != 0:
                        k_values[i] = k_values[j]
                        break
    return k_values

def calculate_gc_k_relation(input_file, output_file, bins_count=250, loess_frac=0.1):
    df = pd.read_csv(input_file, sep='\t', header=None,
                     names=['chr','start','end','gc_content','depth'],
                     dtype={'chr':str,'start':int,'end':int,'gc_content':float,'depth':str})

    df['depth'] = df['depth'].replace('.', '0').fillna('0').astype(float)
    df = df[(df['gc_content'] >= 0) & (df['gc_content'] <= 1)]

    bins = np.linspace(0, 1, bins_count + 1)
    df['bin'] = np.digitize(df['gc_content'], bins) - 1
    df = df[(df['bin'] >= 0) & (df['bin'] < bins_count)]

    bin_stats = df.groupby('bin')['depth'].mean().reindex(range(bins_count), fill_value=0)
    bin_stats.index.name = 'bin'
    bin_stats = bin_stats.reset_index()
    bin_stats.rename(columns={'depth':'mean_depth'}, inplace=True)

    bin_gc_mid = (bins[:-1] + bins[1:]) / 2
    valid_idx = bin_stats['mean_depth'] > 0
    x_fit = bin_gc_mid[valid_idx]
    y_fit = bin_stats.loc[valid_idx, 'mean_depth']

    # Fitting a nonzero bin using LOESS
    loess_result = lowess(y_fit, x_fit, frac=loess_frac, return_sorted=False)

    # Insert the results into the complete loess_fit array
    loess_fit = np.full_like(bin_gc_mid, np.nan)
    loess_fit[valid_idx] = loess_result
    loess_fit = np.clip(loess_fit, 1e-3, None)  # Prevent division by zero

    # Calculate the global average depth using the original depth
    overall_mean = df['depth'].mean()
    print(overall_mean)

    # Calculate the K value within the effective position
    k_loess = np.full(bins_count, np.nan)
    k_loess[valid_idx] = overall_mean / loess_fit[valid_idx]
    k_loess = np.clip(k_loess, None, 2)

    # Fill in the missing K values
    k_loess_corrected = fill_zero_k_values_by_gc_direction(k_loess, bins)

    # Output
    output_df = pd.DataFrame({
        'bin_start': bins[:-1],
        'bin_end': bins[1:],
        'mean_depth': bin_stats['mean_depth'],
        'K_loess': k_loess_corrected
    })

    output_df['bin_start'] = output_df['bin_start'].round(3)
    output_df['bin_end'] = output_df['bin_end'].round(3)
    output_df.to_csv(output_file, sep='\t', index=False)
    print(f"Output saved to {output_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='GC bin depth correction factors with LOESS and GC-based filling.')
    parser.add_argument('input_file', type=str, help='Input file: chr start end gc_content depth (tab-separated)')
    parser.add_argument('output_file', type=str, help='Output file with bin_start, bin_end, mean_depth, K_loess')
    parser.add_argument('--bins_count', type=int, default=250, help='Number of GC bins')
    parser.add_argument('--loess_frac', type=float, default=0.1, help='LOESS smoothing fraction')
    args = parser.parse_args()

    calculate_gc_k_relation(args.input_file, args.output_file, args.bins_count, args.loess_frac)
