# COVO

COVO stands for "Correlations of Vector Orientations" and is a C++ code for computing two-point correlations of 3D vector fields sampled at discrete particle positions.

## Definitions

COVO computes two types of correlations:

$$\eta_{vv}(r) = \langle |{\bf \hat{v}}_{a} \cdot {\bf \hat{v}_b}|^{\gamma}\rangle (r) \ ,$$

and

$$\eta_{vr}(r) = \langle |{\bf \hat{v}}_{a} \cdot {\bf \hat{r}_{ab}}|^{\gamma}\rangle (r) \ ,$$

where ${\bf \hat{v}}_{a}$ and ${\bf \hat{v}}_{b}$ are normalized three-dimensional vectors at positions $a$ and $b$, ${\bf \hat r}_{ab}$ is the normalized distance vector pointing from position $a$ to position $b$, and $\langle \dots \rangle$ denotes the average over all possible vector pairs separated by distance $r = | {\bf r}_{ab} |$. The exponent $\gamma$ is typically set to $1$ or $2$, resulting in $\eta(r)=1/2$ or $\eta(r)=1/3$ for randomly oriented vectors.

## Features

- Both auto- and cross-correlations between two fields can be computed
- Two vectors are assigned to each position, allowing for the computation of multiple correlations in one run
- Operates in both cubical and spherical volumes
- Includes a simple spatial tree structure (using cubical grid cells or HEALPix cells for cubical or spherical volumes respectively) that serves two purposes:
  - Accelerated computation
  - Jackknife sampling for error estimation
- Parallel computation on shared memory using OpenMP


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

## Reference

If you use this code for your research, please provide a link to this repository when publishing the results.

