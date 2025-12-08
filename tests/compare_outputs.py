#!/usr/bin/env python3
"""
Regression test comparison script for covo correlation function outputs.

This script compares two output files from covo, handling:
- Floating point comparisons with tolerance
- Header line skipping
- Column-by-column comparison
- Detailed reporting of differences
"""

import sys
import os
import argparse
import numpy as np

def parse_covo_output(filename):
    """
    Parse a covo output file, skipping header lines and extracting data.
    Returns a dictionary with column names and a numpy array of data.
    """
    data_lines = []
    header = None
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                if line.startswith('#') and header is None:
                    # Extract header (remove '# ' prefix)
                    header = line[2:].split()
                continue
            # Parse data line
            parts = line.split()
            try:
                data_lines.append([float(x) for x in parts])
            except ValueError:
                print(f"Warning: Could not parse line in {filename}: {line}")
                continue
    
    if not data_lines:
        raise ValueError(f"No data found in {filename}")
    
    data_array = np.array(data_lines)
    
    return header, data_array

def compare_files(file1, file2, rtol=1e-5, atol=1e-8, verbose=True):
    """
    Compare two covo output files.
    
    Parameters:
    -----------
    file1 : str
        Path to baseline/reference file
    file2 : str
        Path to test output file
    rtol : float
        Relative tolerance for floating point comparison
    atol : float
        Absolute tolerance for floating point comparison
    verbose : bool
        Print detailed comparison information
    
    Returns:
    --------
    bool : True if files match within tolerance, False otherwise
    """
    if not os.path.exists(file1):
        print(f"Error: Baseline file not found: {file1}")
        return False
    
    if not os.path.exists(file2):
        print(f"Error: Test output file not found: {file2}")
        return False
    
    try:
        header1, data1 = parse_covo_output(file1)
        header2, data2 = parse_covo_output(file2)
    except Exception as e:
        print(f"Error parsing files: {e}")
        return False
    
    # Check if headers match (if present)
    if header1 and header2:
        if header1 != header2:
            print(f"Warning: Headers differ between files")
            if verbose:
                print(f"  Baseline: {header1}")
                print(f"  Test:     {header2}")
    
    # Check dimensions
    if data1.shape != data2.shape:
        print(f"Error: Data dimensions differ")
        print(f"  Baseline: {data1.shape}")
        print(f"  Test:     {data2.shape}")
        return False
    
    # Compare data column by column
    n_rows, n_cols = data1.shape
    all_match = True
    
    for col_idx in range(n_cols):
        col1 = data1[:, col_idx]
        col2 = data2[:, col_idx]
        
        # For integer columns (like bin number, counts), use exact comparison
        if col_idx == 0:  # bin number
            if not np.array_equal(col1, col2):
                print(f"Error: Column {col_idx} (bin) differs")
                diff_mask = col1 != col2
                print(f"  Rows with differences: {np.where(diff_mask)[0]}")
                all_match = False
        elif col_idx == 2:  # counts (unsigned long)
            if not np.array_equal(col1, col2):
                print(f"Error: Column {col_idx} (counts) differs")
                diff_mask = col1 != col2
                diff_rows = np.where(diff_mask)[0]
                print(f"  Rows with differences: {diff_rows[:10]}")  # Show first 10
                if len(diff_rows) > 10:
                    print(f"  ... and {len(diff_rows) - 10} more")
                all_match = False
        else:
            # Floating point columns - use tolerance
            if not np.allclose(col1, col2, rtol=rtol, atol=atol):
                print(f"Warning: Column {col_idx} differs beyond tolerance")
                diff_mask = ~np.isclose(col1, col2, rtol=rtol, atol=atol)
                diff_rows = np.where(diff_mask)[0]
                
                if verbose:
                    print(f"  Rows with differences: {diff_rows[:10]}")
                    if len(diff_rows) > 10:
                        print(f"  ... and {len(diff_rows) - 10} more")
                    
                    # Show some example differences
                    for row_idx in diff_rows[:5]:
                        val1 = col1[row_idx]
                        val2 = col2[row_idx]
                        rel_diff = abs(val1 - val2) / (abs(val1) + atol)
                        print(f"    Row {row_idx}: {val1:.8e} vs {val2:.8e} (rel diff: {rel_diff:.2e})")
                
                all_match = False
    
    if all_match:
        if verbose:
            print(f"✓ Files match within tolerance (rtol={rtol}, atol={atol})")
        return True
    else:
        return False

def main():
    parser = argparse.ArgumentParser(
        description='Compare covo output files for regression testing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s baseline.csv output.csv
  %(prog)s baseline.csv output.csv --rtol 1e-6 --atol 1e-9
  %(prog)s baseline.csv output.csv --quiet
        """
    )
    parser.add_argument('baseline', help='Path to baseline/reference output file')
    parser.add_argument('test_output', help='Path to test output file')
    parser.add_argument('--rtol', type=float, default=1e-5,
                        help='Relative tolerance for floating point comparison (default: 1e-5)')
    parser.add_argument('--atol', type=float, default=1e-8,
                        help='Absolute tolerance for floating point comparison (default: 1e-8)')
    parser.add_argument('--quiet', action='store_true',
                        help='Suppress detailed output (only show pass/fail)')
    
    args = parser.parse_args()
    
    match = compare_files(
        args.baseline,
        args.test_output,
        rtol=args.rtol,
        atol=args.atol,
        verbose=not args.quiet
    )
    
    sys.exit(0 if match else 1)

if __name__ == '__main__':
    main()

