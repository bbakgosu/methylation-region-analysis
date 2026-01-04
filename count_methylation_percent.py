#!/usr/bin/env python3
"""
Regional Methylation Percentage Calculator
특정 genomic region의 전체 methylation percentage를 계산하는 도구
"""

import sys
import argparse
import pandas as pd
from typing import Tuple, Optional
import gzip

def parse_region(region_str: str) -> Tuple[str, int, int]:
    """
    Parse region string in format 'chr:start-end'
    Returns: (chromosome, start, end) with 0-based coordinates
    """
    try:
        chrom, coords = region_str.split(':')
        start, end = coords.split('-')
        # Convert to 0-based if input is 1-based (common in user input)
        # We'll assume user input is 1-based and convert to 0-based
        return chrom, int(start) - 1, int(end)
    except:
        raise ValueError(f"Invalid region format: {region_str}. Use 'chr:start-end'")

def read_bedmethyl(filepath: str) -> pd.DataFrame:
    """
    Read bedMethyl file into a pandas DataFrame
    Handles both regular and gzipped files
    """
    # Define column names based on bedMethyl format
    base_columns = [
        'chrom', 'start', 'end', 'mod_type', 'score', 
        'strand', 'start_dup', 'end_dup', 'color', 'coverage'
    ]
    
    # Check if file is gzipped
    if filepath.endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
    
    # First, check how many columns the file has
    with opener(filepath, mode) as f:
        first_line = f.readline().strip()
        if not first_line or first_line.startswith('#'):
            # Skip header if present
            first_line = f.readline().strip()
        num_columns = len(first_line.split('\t'))
    
    # Add additional columns based on file format
    if num_columns >= 11:
        base_columns.append('percent_modified')
    if num_columns >= 12:
        base_columns.append('n_modified')
    if num_columns >= 13:
        base_columns.append('n_canonical')
    if num_columns >= 14:
        base_columns.append('n_other_mod')
    if num_columns >= 15:
        base_columns.append('n_delete')
    if num_columns >= 16:
        base_columns.append('n_fail')
    if num_columns >= 17:
        base_columns.append('n_diff')
    if num_columns >= 18:
        base_columns.append('n_nocall')
    
    # Read the file
    columns_to_use = base_columns[:num_columns]
    
    df = pd.read_csv(
        filepath,
        sep='\t',
        names=columns_to_use,
        comment='#',
        dtype={
            'chrom': str,
            'start': int,
            'end': int,
            'mod_type': str,
            'coverage': int
        }
    )
    
    # Convert percent_modified to float if it exists
    if 'percent_modified' in df.columns:
        df['percent_modified'] = pd.to_numeric(df['percent_modified'], errors='coerce')
    
    # Convert n_modified and n_canonical to int if they exist
    if 'n_modified' in df.columns:
        df['n_modified'] = pd.to_numeric(df['n_modified'], errors='coerce').fillna(0).astype(int)
    if 'n_canonical' in df.columns:
        df['n_canonical'] = pd.to_numeric(df['n_canonical'], errors='coerce').fillna(0).astype(int)
    
    return df

def calculate_region_methylation(df: pd.DataFrame, chrom: str, start: int, end: int) -> dict:
    """
    Calculate overall methylation percentage for a specific region
    
    Args:
        df: bedMethyl dataframe
        chrom: chromosome name
        start: 0-based start position
        end: 0-based end position (exclusive)
    
    Returns:
        Dictionary with methylation statistics
    """
    # Filter for 5mC modifications only (type 'm')
    df_5mc = df[df['mod_type'] == 'm'].copy()
    
    # Filter for the specified region
    df_region = df_5mc[
        (df_5mc['chrom'] == chrom) & 
        (df_5mc['start'] >= start) & 
        (df_5mc['end'] <= end)
    ].copy()
    
    if len(df_region) == 0:
        return {
            'error': f'No 5mC data found in region {chrom}:{start+1}-{end}',
            'total_positions': 0,
            'total_modified_reads': 0,
            'total_unmodified_reads': 0,
            'total_reads': 0,
            'overall_methylation_percent': 0.0
        }
    
    # Calculate statistics based on available columns
    if 'n_modified' in df_region.columns and 'n_canonical' in df_region.columns:
        # Use actual read counts if available
        total_modified = df_region['n_modified'].sum()
        total_unmodified = df_region['n_canonical'].sum()
        total_reads = total_modified + total_unmodified
        
        # Handle other modifications if column exists
        if 'n_other_mod' in df_region.columns:
            total_other_mod = df_region['n_other_mod'].sum()
        else:
            total_other_mod = 0
            
    elif 'percent_modified' in df_region.columns and 'coverage' in df_region.columns:
        # Calculate from percentage and coverage
        df_region['n_modified_calc'] = (df_region['percent_modified'] / 100.0 * df_region['coverage']).round()
        df_region['n_unmodified_calc'] = df_region['coverage'] - df_region['n_modified_calc']
        
        total_modified = df_region['n_modified_calc'].sum()
        total_unmodified = df_region['n_unmodified_calc'].sum()
        total_reads = df_region['coverage'].sum()
        total_other_mod = 0
    else:
        return {
            'error': 'Insufficient data columns to calculate methylation',
            'total_positions': len(df_region),
            'total_modified_reads': 0,
            'total_unmodified_reads': 0,
            'total_reads': 0,
            'overall_methylation_percent': 0.0
        }
    
    # Calculate overall methylation percentage
    if total_reads > 0:
        overall_methylation_percent = (total_modified / total_reads) * 100
    else:
        overall_methylation_percent = 0.0
    
    # Calculate additional statistics
    positions_covered = len(df_region)
    mean_coverage = df_region['coverage'].mean() if 'coverage' in df_region.columns else 0
    median_coverage = df_region['coverage'].median() if 'coverage' in df_region.columns else 0
    
    # Get strand-specific counts
    plus_strand = df_region[df_region['strand'] == '+']
    minus_strand = df_region[df_region['strand'] == '-']
    
    result = {
        'region': f'{chrom}:{start+1}-{end}',
        'total_positions': positions_covered,
        'total_modified_reads': int(total_modified),
        'total_unmodified_reads': int(total_unmodified),
        'total_reads': int(total_reads),
        'overall_methylation_percent': round(overall_methylation_percent, 2),
        'mean_coverage_per_position': round(mean_coverage, 2),
        'median_coverage_per_position': round(median_coverage, 2),
        'positions_plus_strand': len(plus_strand),
        'positions_minus_strand': len(minus_strand)
    }
    
    if total_other_mod > 0:
        result['total_other_modifications'] = int(total_other_mod)
    
    return result

def print_results(results: dict, verbose: bool = False):
    """
    Print calculation results in a formatted way
    """
    if 'error' in results:
        print(f"Error: {results['error']}", file=sys.stderr)
        return
    
    print("\n" + "="*60)
    print(f"Methylation Analysis for Region: {results['region']}")
    print("="*60)
    
    print(f"\n📊 Overall Statistics:")
    print(f"  • Total CpG positions analyzed: {results['total_positions']:,}")
    print(f"  • Total reads analyzed: {results['total_reads']:,}")
    print(f"  • Mean coverage per position: {results['mean_coverage_per_position']:.1f}")
    print(f"  • Median coverage per position: {results['median_coverage_per_position']:.1f}")
    
    print(f"\n🧬 Methylation Counts:")
    print(f"  • Methylated reads: {results['total_modified_reads']:,}")
    print(f"  • Unmethylated reads: {results['total_unmodified_reads']:,}")
    if 'total_other_modifications' in results:
        print(f"  • Other modifications: {results['total_other_modifications']:,}")
    
    print(f"\n📈 Methylation Percentage:")
    print(f"  • Overall methylation: {results['overall_methylation_percent']:.2f}%")
    
    # Visual representation
    methylation_pct = results['overall_methylation_percent']
    bar_length = 40
    filled = int(bar_length * methylation_pct / 100)
    bar = '█' * filled + '░' * (bar_length - filled)
    print(f"  • Visual: [{bar}] {methylation_pct:.1f}%")
    
    if verbose:
        print(f"\n🔬 Strand Distribution:")
        print(f"  • Plus strand positions: {results['positions_plus_strand']}")
        print(f"  • Minus strand positions: {results['positions_minus_strand']}")
    
    print("\n" + "="*60)

def main():
    parser = argparse.ArgumentParser(
        description='Calculate overall methylation percentage for a genomic region',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with region from command line
  %(prog)s -i sample.bedmethyl -r chr1:14923-15923
  
  # With verbose output
  %(prog)s -i sample.bedmethyl -r chr1:14923-15923 -v
  
  # Gzipped input file
  %(prog)s -i sample.bedmethyl.gz -r chr1:14923-15923
  
  # Output in simple format
  %(prog)s -i sample.bedmethyl -r chr1:14923-15923 --simple
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                      help='Input bedMethyl file (can be gzipped)')
    parser.add_argument('-r', '--region', required=True,
                      help='Genomic region (format: chr:start-end, 1-based coordinates)')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Show detailed output')
    parser.add_argument('--simple', action='store_true',
                      help='Simple output format (just the key numbers)')
    
    args = parser.parse_args()
    
    try:
        # Parse region
        chrom, start, end = parse_region(args.region)
        
        if args.verbose:
            print(f"Loading bedMethyl file: {args.input}")
        
        # Read bedMethyl file
        df = read_bedmethyl(args.input)
        
        if args.verbose:
            print(f"Loaded {len(df):,} total positions")
            print(f"Found {len(df[df['mod_type'] == 'm']):,} 5mC positions")
        
        # Calculate regional methylation
        results = calculate_region_methylation(df, chrom, start, end)
        
        # Print results
        if args.simple:
            # Simple output for easy parsing
            if 'error' not in results:
                print(f"Region: {results['region']}")
                print(f"Positions: {results['total_positions']}")
                print(f"Modified reads: {results['total_modified_reads']}")
                print(f"Unmodified reads: {results['total_unmodified_reads']}")
                print(f"Total reads: {results['total_reads']}")
                print(f"Methylation: {results['overall_methylation_percent']}%")
            else:
                print(f"Error: {results['error']}", file=sys.stderr)
        else:
            print_results(results, verbose=args.verbose)
            
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()