# covo

computes the average inner product of vector pairs as a function of their distance

---

## compile
$ cd src
$ make


## run

option 0:   all parameters are specified by the parameter file <br />
./covo covo.params <br />

option 1:   catalogue_1 is input argument, catalogue_2 is copied from catalogue_1 <br />
./covo covo.params catalogue_1.csv <br />

option 2:   catalogue_1 and catalogue_2 are input arguments <br />
./covo covo.params catalogue_1.csv catalogue_2.csv <br />

option 3:   suffix for output file as last argument <br />
./covo covo.params catalogue_1.csv catalogue_2.csv suffix_for_output_file <br />


## parameters
Parameters are set and described in covo.params


## input catalogues
Code searches for pairs of objects belonging to two different input catalogues. <br />

Input files should be ASCII and contain information on the position and <br />
components of two vectors for each object, i.e. <br />
x, y, z, vec_a_x, vec_a_y, vec_a_z, vec_b_x, vec_b_y, vec_b_z <br />

Both input catalogues must have the same format. <br />
Columns and column delimiter are set in parameter file.

