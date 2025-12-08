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
BASELINE_DIR="tests/baseline"
OUTPUT_DIR="tests/output"
COMPARE_SCRIPT="tests/compare_outputs.py"

# Test cases: (test_name, param_file, cat1, cat2, suffix)
# suffix is used as the 4th argument to covo and in output filename
# Note: cat1 and cat2 are also specified in param files, but can be overridden via command line
declare -a TEST_CASES=(
    "box_mode:tests/covo.params.test_box_mode:./tests/data/box_mode_cat1.csv:./tests/data/box_mode_cat2.csv:regression_test_box"
    "shell_mode:tests/covo.params.test_shell_mode:./tests/data/shell_mode_cat1.csv:./tests/data/shell_mode_cat2.csv:regression_test_shell"
)

echo "=========================================="
echo "COVO Regression Test Suite"
echo "=========================================="
echo ""

# Check if covo executable exists
if [ ! -f "./covo" ]; then
    echo -e "${RED}Error: covo executable not found. Please build first with 'make'${NC}"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$BASELINE_DIR"

# Track test results
FAILED_TESTS=()
PASSED_TESTS=()

# Function to run a single test case
run_test_case() {
    local test_name="$1"
    local param_file="$2"
    local cat1="$3"
    local cat2="$4"
    local suffix="$5"
    
    # Construct output filename (matches auto mode: dir_out + prefix + "_" + suffix + "." + extension)
    local test_output="${OUTPUT_DIR}/w_${suffix}.csv"
    local baseline_output="${BASELINE_DIR}/w_${suffix}.csv"
    
    echo "----------------------------------------"
    echo "Test: $test_name"
    echo "----------------------------------------"
    echo ""
    
    # Check if parameter file exists
    if [ ! -f "$param_file" ]; then
        echo -e "${RED}Error: Parameter file not found: $param_file${NC}"
        FAILED_TESTS+=("$test_name (missing param file)")
        return 1
    fi
    
    # Check if baseline exists
    if [ ! -f "$baseline_output" ]; then
        echo -e "${YELLOW}Warning: Baseline file not found: $baseline_output${NC}"
        echo -e "${YELLOW}This is expected for the first run.${NC}"
        echo ""
        echo "Running test to generate baseline..."
        echo ""
        
        # Run the test
        ./covo "$param_file" "$cat1" "$cat2" "$suffix"
        
        # Copy to baseline
        cp "$test_output" "$baseline_output"
        
        echo ""
        echo -e "${GREEN}✓ Baseline created: $baseline_output${NC}"
        echo -e "${YELLOW}Please review the baseline and commit it if correct.${NC}"
        echo ""
        PASSED_TESTS+=("$test_name (baseline generated)")
        return 0
    fi
    
    # Run the test
    echo "Running covo..."
    echo "  Parameter file: $param_file"
    echo "  Input: $cat1, $cat2"
    echo "  Output: $test_output"
    echo ""
    
    ./covo "$param_file" "$cat1" "$cat2" "$suffix"
    
    # Compare outputs
    echo ""
    echo "Comparing outputs..."
    echo "  Baseline: $baseline_output"
    echo "  Test:     $test_output"
    echo ""
    
    if python3 "$COMPARE_SCRIPT" "$baseline_output" "$test_output"; then
        echo ""
        echo -e "${GREEN}✓ Test PASSED: $test_name${NC}"
        echo ""
        PASSED_TESTS+=("$test_name")
        return 0
    else
        echo ""
        echo -e "${RED}✗ Test FAILED: $test_name${NC}"
        echo ""
        echo "If the differences are expected (e.g., due to intentional code changes),"
        echo "update the baseline by running:"
        echo "  cp $test_output $baseline_output"
        echo ""
        FAILED_TESTS+=("$test_name")
        return 1
    fi
}

# Run all test cases
for test_case in "${TEST_CASES[@]}"; do
    IFS=':' read -r test_name param_file cat1 cat2 suffix <<< "$test_case"
    run_test_case "$test_name" "$param_file" "$cat1" "$cat2" "$suffix" || true
done

# Print summary
echo "=========================================="
echo "Test Summary"
echo "=========================================="
echo ""

if [ ${#PASSED_TESTS[@]} -gt 0 ]; then
    echo -e "${GREEN}Passed tests (${#PASSED_TESTS[@]}):${NC}"
    for test in "${PASSED_TESTS[@]}"; do
        echo "  ✓ $test"
    done
    echo ""
fi

if [ ${#FAILED_TESTS[@]} -gt 0 ]; then
    echo -e "${RED}Failed tests (${#FAILED_TESTS[@]}):${NC}"
    for test in "${FAILED_TESTS[@]}"; do
        echo "  ✗ $test"
    done
    echo ""
    exit 1
fi

if [ ${#PASSED_TESTS[@]} -eq 0 ] && [ ${#FAILED_TESTS[@]} -eq 0 ]; then
    echo -e "${YELLOW}No tests were run (all baselines were generated)${NC}"
    echo ""
    exit 0
fi

# Generate comparison plot if both test outputs exist
BOX_BASELINE="${BASELINE_DIR}/w_regression_test_box.csv"
BOX_TEST="${OUTPUT_DIR}/w_regression_test_box.csv"
SHELL_BASELINE="${BASELINE_DIR}/w_regression_test_shell.csv"
SHELL_TEST="${OUTPUT_DIR}/w_regression_test_shell.csv"

if [ -f "$BOX_BASELINE" ] && [ -f "$BOX_TEST" ] && \
   [ -f "$SHELL_BASELINE" ] && [ -f "$SHELL_TEST" ]; then
    echo ""
    echo "Generating comparison plot..."
    if python3 "$COMPARE_SCRIPT" plot \
        --baseline-box "$BOX_BASELINE" \
        --test-box "$BOX_TEST" \
        --baseline-shell "$SHELL_BASELINE" \
        --test-shell "$SHELL_TEST" \
        --output "${OUTPUT_DIR}/comparison_plot.png"; then
        echo -e "${GREEN}✓ Comparison plot saved to ${OUTPUT_DIR}/comparison_plot.png${NC}"
    else
        echo -e "${YELLOW}Warning: Could not generate comparison plot${NC}"
    fi
fi

echo -e "${GREEN}=========================================="
echo "✓ All regression tests PASSED"
echo "==========================================${NC}"
exit 0
