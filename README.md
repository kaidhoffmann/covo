# COVO

COVO (Correlations of Vector Orientations) is a C++ code for computing correlation functions of 3D vector fields sampled at discrete particle positions.

## Features

- Measures correlations of discrete vector fields as the average inner product:
  - between vectors at two different positions (vector-vector correlations)
  - between a vector at one position and the separation vector pointing to another position (vector-position correlations)
- Correlations are computed as a function of the distance between two positions
- Both auto- and cross-correlations between two fields can be computed
- Two vectors are assigned to each position, allowing for the computation of multiple correlations in one run
- Operates in both cubical and spherical volumes
- Includes a tree structure that serves two purposes:
  - Accelerated computation
  - Jackknife sampling for error estimation, using cubical grid cells or HEALPix cells for cubical or spherical volumes respectively
- Parallel computation using OpenMP

## Citation

If you use this code for your research, please provide a link to this repository when publishing the results.

---

## Compilation

```bash
make -C src
```

The executable `covo` will be created in the project root directory.

## Usage

### Basic Usage

COVO searches for pairs of objects belonging to two different input catalogues. The catalogues can be specified in the parameter file or provided as command-line arguments.

**Option 1:** All parameters are specified by the parameter file
```bash
./covo covo.params
```

**Option 2:** Catalogue_1 is provided as input argument, catalogue_2 is copied from catalogue_1
```bash
./covo covo.params catalogue_1.csv
```

**Option 3:** Both catalogues are provided as input arguments
```bash
./covo covo.params catalogue_1.csv catalogue_2.csv
```

**Option 4:** Output file or suffix for output file as last argument (other options are set in parameter file)
```bash
./covo covo.params catalogue_1.csv catalogue_2.csv output_file.csv
```

## Parameters

Parameters are set and described in `covo.params`. See the parameter file for detailed documentation of all available options.

## File Format

**Input files** are CSV files containing information on the position and components of two vectors for each object:

```
x, y, z, vec_a_x, vec_a_y, vec_a_z, vec_b_x, vec_b_y, vec_b_z
```

Both input catalogues must have the same format. Column numbers and delimiters are configured in the parameter file.

**Output files** are CSV files with one row per distance bin. The output format (header presence, delimiter, column selection) can be configured via parameters in the parameter file (`header_out`, `delim_out`).

## Regression Testing

Tools for regression testing are available to ensure code changes do not alter numerical results. See [`tests/README.md`](tests/README.md) for details.

Run regression tests with:
```bash
make -C src test
```
