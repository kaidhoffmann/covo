# Regression Testing for COVO

This directory contains the regression testing framework for the COVO correlation function computation code.

## Overview

Regression tests ensure that code changes do not alter the numerical results of correlation function computations. The framework:

1. Runs COVO with predefined test inputs
2. Compares outputs against baseline reference files
3. Reports differences with detailed diagnostics

## Directory Structure

```
tests/
├── README.md                    # This file
├── covo.params.test_box_mode   # Box mode test parameter file
├── covo.params.test_shell_mode # Shell mode test parameter file
├── compare_outputs.py          # Python script for comparing outputs
├── run_regression_tests.sh     # Main test runner script
├── data/                        # Test input catalogues
│   ├── box_mode_cat1.csv       # Box mode catalogue 1
│   ├── box_mode_cat2.csv       # Box mode catalogue 2
│   ├── shell_mode_cat1.csv     # Shell mode catalogue 1
│   └── shell_mode_cat2.csv     # Shell mode catalogue 2
├── baseline/                    # Reference output files (committed to git)
│   └── regression_test.csv     # Baseline output for main test
└── output/                      # Test output files (git-ignored)
    └── regression_test.csv     # Current test run output
```

## Quick Start

### Running Tests

From the project root directory:

```bash
# Run regression tests (builds automatically if needed)
make -C src test
# OR
bash tests/run_regression_tests.sh
```

### First Run (Creating Baseline)

On the first run, the baseline file doesn't exist. The test script will:

1. Run the test
2. Automatically create the baseline from the output
3. Prompt you to review and commit it

**Important:** Review the generated baseline carefully before committing it to ensure it represents the correct expected output.

### Updating Baseline

If you make intentional changes that affect the output (e.g., algorithm improvements, bug fixes), update the baseline:

```bash
# Run tests (they will fail)
make -C src test

# If the differences are expected, update the baseline
cp tests/output/regression_test.csv tests/baseline/regression_test.csv

# Verify the update
make -C src test
```

## Customizing Tests

To add more test cases:

1. Create additional parameter files (e.g., `covo.params.test_new`)
2. Add test cases to `run_regression_tests.sh`
3. Generate baseline outputs for each test case

## Output Comparison

The `compare_outputs.py` script compares outputs with:

- **Exact matching** for integer columns (bin numbers, counts)
- **Tolerance-based comparison** for floating-point columns (correlation values, errors)
  - Default relative tolerance: 1e-5
  - Default absolute tolerance: 1e-8

### Adjusting Tolerance

If you need different tolerances:

```bash
python3 tests/compare_outputs.py \
    tests/baseline/regression_test.csv \
    tests/output/regression_test.csv \
    --rtol 1e-6 \
    --atol 1e-9
```

## Integration with Development Workflow

### Before Committing Code

Always run regression tests before committing:

```bash
make -C src test
```

### CI/CD Integration

The test script exits with:
- **0** on success
- **1** on failure

This makes it suitable for CI/CD pipelines. See [`.github/workflows/regression.yml`](../.github/workflows/regression.yml) for the GitHub Actions workflow.

## Adding New Test Cases

To add a new test case:

1. **Create test parameter file**:
   ```bash
   cp tests/covo.params.test_box_mode tests/covo.params.test_new
   # Edit parameters as needed
   ```

2. **Add to test runner**:
   Edit `run_regression_tests.sh` to include the new test case

3. **Generate baseline**:
   Run the test once to generate the baseline, then move it to `baseline/`

4. **Document**:
   Add a comment in the test runner explaining what the test validates
