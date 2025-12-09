# COVO

COVO (Correlations of Vector Orientations) is a C++ code for computing correlation functions of 3D vector fields sampled at discrete particle positions.

## Features

- Measures correlations of discrete vector fields as the average inner product:
  - between vectors at two different positions (vector-vector correlations)
  - between a vector at one position and the separation vector pointing to another position (vector-position correlations)
- Correlations are computed as a function of the distance between two positions
- Both auto- and cross-correlations between two fields can be computed
- Two vectors are assigned to each position, allowing for the computation of multiple correlations in one run for increased efficiency
- Operates in both cubical and spherical volumes
- Includes a tree structure that serves two purposes:
  - Accelerated computation
  - Jackknife sampling for error estimation, using cubical grid cells or HEALPix cells for cubical or spherical volumes respectively

## Citation

If you use this code for your research, please provide a link to this repository when publishing the results.

---

## Compilation

```bash
cd src
make
```

The executable `covo` will be created in the project root directory.

## Usage

### Basic Usage

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

## Input Catalogues

The code searches for pairs of objects belonging to two different input catalogues.

### Format

Input files should be ASCII and contain information on the position and components of two vectors for each object:

```
x, y, z, vec_a_x, vec_a_y, vec_a_z, vec_b_x, vec_b_y, vec_b_z
```

Both input catalogues must have the same format. Columns and column delimiter are set in the parameter file.

## Regression Testing

Tools for regression testing are available to ensure code changes do not alter numerical results. See [`REGRESSION_TESTING.md`](REGRESSION_TESTING.md) for details.

Run regression tests with:
```bash
make -C src test
```
