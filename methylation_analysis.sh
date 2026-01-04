#!/bin/bash

# ============================================================
# Methylation Analysis Pipeline
# Analyzes methylation in a specific genomic region from BAM file
# ============================================================

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored messages
print_message() {
    local color=$1
    local message=$2
    echo -e "${color}${message}${NC}"
}

# Function to print usage
usage() {
    cat << EOF
Usage: $(basename "$0") -b <bam_file> -r <region> [options]

Methylation Analysis Pipeline
Analyzes 5mC methylation in a specific genomic region from BAM files

Required Arguments:
    -b, --bam FILE        Input BAM file with methylation tags
    -r, --region REGION   Genomic region (format: chr:start-end, e.g., chr1:14923-15923)
                          Commas allowed for readability (e.g., chr1:1,000,000-2,000,000)

Optional Arguments:
    -o, --output DIR      Output directory (default: ./output)
    -t, --threads NUM     Number of threads for modkit (default: 8)
    -s, --script PATH     Path to count_methylation_percent.py (default: ./count_methylation_percent.py)
    -m, --modkit PATH     Path to modkit executable (default: modkit)
    --simple              Use simple output format for methylation stats
    -v, --verbose         Enable verbose output
    -h, --help            Show this help message and exit

Examples:
    # Basic usage
    $(basename "$0") -b sample.bam -r chr1:14923-15923

    # With commas in region
    $(basename "$0") -b sample.bam -r "chr1:1,000,000-2,000,000"

    # With custom output directory and threads
    $(basename "$0") -b sample.bam -r chr1:14923-15923 -o results -t 16

    # With custom script path
    $(basename "$0") -b sample.bam -r chr1:14923-15923 -s /path/to/count_methylation_percent.py

    # Verbose mode with simple output
    $(basename "$0") -b sample.bam -r chr1:14923-15923 -v --simple

EOF
    exit 1
}

# Function to check if command exists
check_command() {
    local cmd=$1
    if ! command -v "$cmd" &> /dev/null; then
        print_message "$RED" "Error: $cmd is not installed or not in PATH"
        exit 1
    fi
}

# Function to check if file exists
check_file() {
    local file=$1
    local description=$2
    if [[ ! -f "$file" ]]; then
        print_message "$RED" "Error: $description not found: $file"
        exit 1
    fi
}

# Function to validate region format
validate_region() {
    local region=$1
    # Allow commas in numbers for readability (e.g., chr1:1,000,000-2,000,000)
    if [[ ! $region =~ ^[a-zA-Z0-9]+:[0-9,]+-[0-9,]+$ ]]; then
        print_message "$RED" "Error: Invalid region format. Use chr:start-end (e.g., chr1:14923-15923 or chr1:1,000,000-2,000,000)"
        exit 1
    fi
}

# Function to clean region format (remove commas for processing)
clean_region() {
    local region=$1
    # Remove commas from numbers
    echo "$region" | tr -d ','
}

# Default values
BAM_FILE=""
REGION=""
OUTPUT_DIR="./output"
THREADS=8
PYTHON_SCRIPT="./count_methylation_percent.py"
MODKIT_CMD="modkit"
SIMPLE_OUTPUT=""
VERBOSE=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -b|--bam)
            BAM_FILE="$2"
            shift 2
            ;;
        -r|--region)
            REGION="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -s|--script)
            PYTHON_SCRIPT="$2"
            shift 2
            ;;
        -m|--modkit)
            MODKIT_CMD="$2"
            shift 2
            ;;
        --simple)
            SIMPLE_OUTPUT="--simple"
            shift
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            print_message "$RED" "Error: Unknown option $1"
            usage
            ;;
    esac
done

# Validate required arguments
if [[ -z "$BAM_FILE" ]] || [[ -z "$REGION" ]]; then
    print_message "$RED" "Error: Missing required arguments"
    usage
fi

# Validate inputs
check_file "$BAM_FILE" "BAM file"
check_file "$PYTHON_SCRIPT" "Python analysis script (count_methylation_percent.py)"
validate_region "$REGION"

# Clean region format (remove commas if present)
REGION_CLEAN=$(clean_region "$REGION")
REGION_DISPLAY="$REGION"  # Keep original for display

# Check for required commands
check_command "$MODKIT_CMD"
check_command "python3"

# Extract sample name from BAM file
SAMPLE_NAME=$(basename "$BAM_FILE" .bam)
SAMPLE_NAME=$(basename "$SAMPLE_NAME" .sorted)
SAMPLE_NAME=$(basename "$SAMPLE_NAME" .aligned)

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Set up file paths
BEDMETHYL_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_bedmethyl.bed"
LOG_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_methylation_analysis.log"
SUMMARY_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_methylation_summary.txt"

# Print header
print_message "$BLUE" "\n=========================================="
print_message "$BLUE" "   Methylation Analysis Pipeline"
print_message "$BLUE" "=========================================="
print_message "$GREEN" "Sample: $SAMPLE_NAME"
print_message "$GREEN" "Region: $REGION_DISPLAY"
if [[ "$REGION_DISPLAY" != "$REGION_CLEAN" ]]; then
    print_message "$YELLOW" "Processing as: $REGION_CLEAN"
fi
print_message "$GREEN" "Output: $OUTPUT_DIR"
print_message "$BLUE" "=========================================="

# Start logging
{
    echo "Methylation Analysis Pipeline Log"
    echo "Date: $(date)"
    echo "Sample: $SAMPLE_NAME"
    echo "BAM file: $BAM_FILE"
    echo "Region: $REGION_DISPLAY (processed as: $REGION_CLEAN)"
    echo "Threads: $THREADS"
    echo "Python script: $PYTHON_SCRIPT"
    echo "==========================================\n"
} > "$LOG_FILE"

# Step 1: Run modkit pileup
print_message "$YELLOW" "\n[Step 1/2] Running modkit pileup..."

# Build modkit command (use cleaned region without commas)
MODKIT_CMD_FULL="$MODKIT_CMD pileup -t $THREADS --region $REGION_CLEAN $BAM_FILE $BEDMETHYL_FILE"

if $VERBOSE; then
    echo "Command: $MODKIT_CMD_FULL"
fi

# Run modkit
if $MODKIT_CMD_FULL 2>> "$LOG_FILE"; then
    print_message "$GREEN" "✓ modkit pileup completed successfully"
    echo "modkit pileup: SUCCESS" >> "$LOG_FILE"
    
    # Check if bedMethyl file was created and has content
    if [[ ! -s "$BEDMETHYL_FILE" ]]; then
        print_message "$YELLOW" "Warning: bedMethyl file is empty. Region might not contain methylation data."
        echo "Warning: Empty bedMethyl file" >> "$LOG_FILE"
    else
        LINE_COUNT=$(wc -l < "$BEDMETHYL_FILE")
        print_message "$GREEN" "  Generated bedMethyl file with $LINE_COUNT lines"
        echo "bedMethyl lines: $LINE_COUNT" >> "$LOG_FILE"
        
        if $VERBOSE; then
            # Show first few lines of bedMethyl file
            print_message "$BLUE" "\n  First 5 lines of bedMethyl file:"
            head -n 5 "$BEDMETHYL_FILE" | while IFS= read -r line; do
                echo "    $line"
            done
        fi
    fi
else
    print_message "$RED" "✗ modkit pileup failed. Check $LOG_FILE for details"
    echo "modkit pileup: FAILED" >> "$LOG_FILE"
    exit 1
fi

# Step 2: Analyze methylation using Python script
print_message "$YELLOW" "\n[Step 2/2] Analyzing regional methylation..."

# Build Python command (use cleaned region without commas)
PYTHON_CMD="python3 $PYTHON_SCRIPT -i $BEDMETHYL_FILE -r $REGION_CLEAN $SIMPLE_OUTPUT"

if $VERBOSE; then
    echo "Command: $PYTHON_CMD"
fi

# Run Python analysis and capture output
ANALYSIS_OUTPUT=$(eval "$PYTHON_CMD" 2>&1)
ANALYSIS_EXIT_CODE=$?

if [[ $ANALYSIS_EXIT_CODE -eq 0 ]]; then
    print_message "$GREEN" "✓ Methylation analysis completed"
    
    # Display analysis output
    echo -e "$ANALYSIS_OUTPUT"
    
    # Save analysis to log
    {
        echo -e "\n==========================================\n"
        echo "Methylation Analysis Results:"
        echo -e "$ANALYSIS_OUTPUT"
    } >> "$LOG_FILE"
    
    # Save summary to separate file
    {
        echo "Sample: $SAMPLE_NAME"
        echo "BAM file: $BAM_FILE"
        echo "Region: $REGION_DISPLAY"
        echo "Analysis Date: $(date)"
        echo "modkit version: $($MODKIT_CMD --version 2>/dev/null || echo 'unknown')"
        echo ""
        echo "==========================================\n"
        echo "$ANALYSIS_OUTPUT"
    } > "$SUMMARY_FILE"
    
    print_message "$GREEN" "\n✓ Results saved to:"
    print_message "$GREEN" "  • bedMethyl file: $BEDMETHYL_FILE"
    print_message "$GREEN" "  • Summary: $SUMMARY_FILE"
    print_message "$GREEN" "  • Log: $LOG_FILE"
    
    # If verbose, show file sizes
    if $VERBOSE; then
        print_message "$BLUE" "\nFile sizes:"
        ls -lh "$BEDMETHYL_FILE" "$SUMMARY_FILE" "$LOG_FILE" | while IFS= read -r line; do
            echo "  $line"
        done
    fi
else
    print_message "$RED" "✗ Methylation analysis failed"
    echo -e "$ANALYSIS_OUTPUT"
    echo "Python analysis: FAILED" >> "$LOG_FILE"
    echo "$ANALYSIS_OUTPUT" >> "$LOG_FILE"
    exit 1
fi

# Final summary in log
{
    echo -e "\n==========================================\n"
    echo "Pipeline completed successfully at $(date)"
    echo "Total execution time: ${SECONDS} seconds"
} >> "$LOG_FILE"

print_message "$GREEN" "\n=========================================="
print_message "$GREEN" "   Pipeline completed successfully!"
print_message "$GREEN" "   Execution time: ${SECONDS} seconds"
print_message "$GREEN" "=========================================="

# Optional: Print quick summary
if [[ -z "$SIMPLE_OUTPUT" ]]; then
    # Extract key metrics from output if not in simple mode
    METHYLATION_PERCENT=$(echo "$ANALYSIS_OUTPUT" | grep -oP "Overall methylation: \K[0-9.]+(?=%)" | head -1)
    TOTAL_READS=$(echo "$ANALYSIS_OUTPUT" | grep -oP "Total reads analyzed: \K[0-9,]+" | head -1 | tr -d ',')
    
    if [[ -n "$METHYLATION_PERCENT" ]]; then
        print_message "$BLUE" "\n📊 Quick Summary:"
        print_message "$BLUE" "  • Region: $REGION_DISPLAY"
        print_message "$BLUE" "  • Methylation: ${METHYLATION_PERCENT}%"
        if [[ -n "$TOTAL_READS" ]]; then
            print_message "$BLUE" "  • Total reads: $(printf "%'d" $TOTAL_READS)"
        fi
    fi
fi