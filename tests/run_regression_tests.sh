#!/bin/bash
# Regression test runner for covo
# This script runs the regression tests and compares outputs against baseline

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"

cd "$PROJECT_ROOT"

# Configuration
TEST_PARAMS="tests/covo.params.test"
BASELINE_DIR="tests/baseline"
OUTPUT_DIR="tests/output"
COMPARE_SCRIPT="tests/compare_outputs.py"

# Test output filename (matches what's set in covo.params.test)
# With fname_out=auto, output is: dir_out + prefix + "_" + suffix + "." + extension
# = ./tests/output/ + w + _ + regression_test + . + csv
TEST_OUTPUT="${OUTPUT_DIR}/w_regression_test.csv"
BASELINE_OUTPUT="${BASELINE_DIR}/w_regression_test.csv"

echo "=========================================="
echo "COVO Regression Test Suite"
echo "=========================================="
echo ""

# Check if covo executable exists
if [ ! -f "./covo" ]; then
    echo -e "${RED}Error: covo executable not found. Please build first with 'make'${NC}"
    exit 1
fi

# Check if test parameter file exists
if [ ! -f "$TEST_PARAMS" ]; then
    echo -e "${RED}Error: Test parameter file not found: $TEST_PARAMS${NC}"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if baseline exists
if [ ! -f "$BASELINE_OUTPUT" ]; then
    echo -e "${YELLOW}Warning: Baseline file not found: $BASELINE_OUTPUT${NC}"
    echo -e "${YELLOW}This is expected for the first run.${NC}"
    echo ""
    echo "Running test to generate baseline..."
    echo ""
    
    # Run the test
    ./covo "$TEST_PARAMS" \
        ./catalogues/random_1.csv \
        ./catalogues/random_2.csv \
        regression_test
    
    # Copy to baseline
    mkdir -p "$BASELINE_DIR"
    cp "$TEST_OUTPUT" "$BASELINE_OUTPUT"
    
    echo ""
    echo -e "${GREEN}✓ Baseline created: $BASELINE_OUTPUT${NC}"
    echo -e "${YELLOW}Please review the baseline and commit it if correct.${NC}"
    exit 0
fi

# Run the test
echo "Running covo with test parameters..."
echo "  Input: ./catalogues/random_1.csv, ./catalogues/random_2.csv"
echo "  Output: $TEST_OUTPUT"
echo ""

./covo "$TEST_PARAMS" \
    ./catalogues/random_1.csv \
    ./catalogues/random_2.csv \
    regression_test

# Compare outputs
echo ""
echo "Comparing outputs..."
echo "  Baseline: $BASELINE_OUTPUT"
echo "  Test:     $TEST_OUTPUT"
echo ""

if python3 "$COMPARE_SCRIPT" "$BASELINE_OUTPUT" "$TEST_OUTPUT"; then
    echo ""
    echo -e "${GREEN}=========================================="
    echo "✓ All regression tests PASSED"
    echo "==========================================${NC}"
    exit 0
else
    echo ""
    echo -e "${RED}=========================================="
    echo "✗ Regression tests FAILED"
    echo "==========================================${NC}"
    echo ""
    echo "If the differences are expected (e.g., due to intentional code changes),"
    echo "update the baseline by running:"
    echo "  cp $TEST_OUTPUT $BASELINE_OUTPUT"
    exit 1
fi

