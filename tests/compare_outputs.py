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
try:
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

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

def plot_comparison(baseline_box, test_box, baseline_shell, test_shell, output_file='tests/output/comparison_plot.png', x_shift=0.02):
    """
    Generate a comparison plot with two panels (box mode and shell mode).
    
    Parameters:
    -----------
    baseline_box : str
        Path to box mode baseline file
    test_box : str
        Path to box mode test output file
    baseline_shell : str
        Path to shell mode baseline file
    test_shell : str
        Path to shell mode test output file
    output_file : str
        Path to save the plot
    x_shift : float
        Shift factor for x-values to separate baseline and test points
        Baseline: x * (1 - x_shift), Test: x * (1 + x_shift)
    """
    if not HAS_MATPLOTLIB:
        print("Warning: matplotlib not available, skipping plot generation")
        return
    
    # Parse all files
    try:
        _, data_box_base = parse_covo_output(baseline_box)
        _, data_box_test = parse_covo_output(test_box)
        _, data_shell_base = parse_covo_output(baseline_shell)
        _, data_shell_test = parse_covo_output(test_shell)
    except Exception as e:
        print(f"Error parsing files for plotting: {e}")
        return
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot box mode (left panel)
    r_box_base = data_box_base[:, 1]  # Column 1: r
    r12_v1a_box_base = data_box_base[:, 3]  # Column 3: r12_v1a
    r12_v1a_std_box_base = data_box_base[:, 4]  # Column 4: r12_v1a_std
    r12_v1b_box_base = data_box_base[:, 5]  # Column 5: r12_v1b
    r12_v1b_std_box_base = data_box_base[:, 6]  # Column 6: r12_v1b_std
    
    r_box_test = data_box_test[:, 1]
    r12_v1a_box_test = data_box_test[:, 3]
    r12_v1a_std_box_test = data_box_test[:, 4]
    r12_v1b_box_test = data_box_test[:, 5]
    r12_v1b_std_box_test = data_box_test[:, 6]
    
    # Baseline plots (solid lines) - offset x by (1 - x_shift)
    ax1.errorbar(r_box_base * (1 - x_shift), r12_v1a_box_base, yerr=r12_v1a_std_box_base,
                 fmt='o-', color='red', label='baseline r12_v1a', capsize=3)
    ax1.errorbar(r_box_base * (1 - x_shift), r12_v1b_box_base, yerr=r12_v1b_std_box_base,
                 fmt='o-', color='red', label='baseline r12_v1b', capsize=3)
    
    # Test plots (x markers) - offset x by (1 + x_shift)
    ax1.errorbar(r_box_test * (1 + x_shift), r12_v1a_box_test, yerr=r12_v1a_std_box_test,
                 fmt='x', color='blue', label='test r12_v1a', capsize=2)
    ax1.errorbar(r_box_test * (1 + x_shift), r12_v1b_box_test, yerr=r12_v1b_std_box_test,
                 fmt='x', color='blue', label='test r12_v1b', capsize=2)
    
    # Reference line at 0.5
    ax1.axhline(y=0.5, color='gray', linestyle=':', alpha=0.5, linewidth=1)
    
    ax1.set_xlabel('r', fontsize=12)
    ax1.set_ylabel('Correlation', fontsize=12)
    ax1.set_title('Box Mode', fontsize=14, fontweight='bold')
    ax1.set_xscale('log')
    ax1.legend(fontsize=9, loc='best')
    ax1.grid(True, alpha=0.3)
    
    # Plot shell mode (right panel)
    r_shell_base = data_shell_base[:, 1]
    r12_v1a_shell_base = data_shell_base[:, 3]
    r12_v1a_std_shell_base = data_shell_base[:, 4]
    r12_v1b_shell_base = data_shell_base[:, 5]
    r12_v1b_std_shell_base = data_shell_base[:, 6]
    
    r_shell_test = data_shell_test[:, 1]
    r12_v1a_shell_test = data_shell_test[:, 3]
    r12_v1a_std_shell_test = data_shell_test[:, 4]
    r12_v1b_shell_test = data_shell_test[:, 5]
    r12_v1b_std_shell_test = data_shell_test[:, 6]
    
    # Baseline plots (solid lines) - offset x by (1 - x_shift)
    ax2.errorbar(r_shell_base * (1 - x_shift), r12_v1a_shell_base, yerr=r12_v1a_std_shell_base,
                 fmt='o-', color='red', label='baseline r12_v1a', capsize=3)
    ax2.errorbar(r_shell_base * (1 - x_shift), r12_v1b_shell_base, yerr=r12_v1b_std_shell_base,
                 fmt='o-', color='red', label='baseline r12_v1b', capsize=3)
    
    # Test plots (x markers) - offset x by (1 + x_shift)
    ax2.errorbar(r_shell_test * (1 + x_shift), r12_v1a_shell_test, yerr=r12_v1a_std_shell_test,
                 fmt='x', color='blue', label='test r12_v1a', capsize=2)
    ax2.errorbar(r_shell_test * (1 + x_shift), r12_v1b_shell_test, yerr=r12_v1b_std_shell_test,
                 fmt='x', color='blue', label='test r12_v1b', capsize=2)
    
    # Reference line at 0.5
    ax2.axhline(y=0.5, color='gray', linestyle=':', alpha=0.5, linewidth=1)
    
    ax2.set_xlabel('r', fontsize=12)
    ax2.set_ylabel('Correlation', fontsize=12)
    ax2.set_title('Shell Mode', fontsize=14, fontweight='bold')
    ax2.set_xscale('log')
    ax2.legend(fontsize=9, loc='best')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else '.', exist_ok=True)
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Plot saved to: {output_file}")
    plt.close()

if __name__ == '__main__':
    # Check if we're being called with plot command
    if len(sys.argv) > 1 and sys.argv[1] == 'plot':
        # Remove 'plot' from argv before parsing
        sys.argv.pop(1)
        # Plot mode: compare all baseline vs test files
        parser = argparse.ArgumentParser(description='Generate comparison plot for regression tests')
        parser.add_argument('--baseline-box', default='tests/baseline/w_regression_test_box.csv',
                          help='Path to box mode baseline file')
        parser.add_argument('--test-box', default='tests/output/w_regression_test_box.csv',
                          help='Path to box mode test output file')
        parser.add_argument('--baseline-shell', default='tests/baseline/w_regression_test_shell.csv',
                          help='Path to shell mode baseline file')
        parser.add_argument('--test-shell', default='tests/output/w_regression_test_shell.csv',
                          help='Path to shell mode test output file')
        parser.add_argument('--output', default='tests/output/comparison_plot.png',
                          help='Output file for the plot')
        parser.add_argument('--x-shift', type=float, default=0.02,
                          help='Shift factor for x-values to separate baseline and test points (default: 0.05)')
        
        args = parser.parse_args()
        plot_comparison(args.baseline_box, args.test_box, args.baseline_shell, args.test_shell, args.output, args.x_shift)
    else:
        # Normal comparison mode
        main()

