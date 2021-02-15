# covo

computes the average inner product of vector pairs as a function of their distance

---

## compile
$ cd src  
$ make


## run

option 0:   all parameters are specified by the parameter file  
./covo covo.params  

option 1:   catalogue_1 is input argument, catalogue_2 is copied from catalogue_1  
./covo covo.params catalogue_1.csv  

option 2:   catalogue_1 and catalogue_2 are input arguments  
./covo covo.params catalogue_1.csv catalogue_2.csv  

option 3:   output file or suffix for output file as last argument (options can are set in parameter file)  
./covo covo.params catalogue_1.csv catalogue_2.csv output_file.csv


## parameters
Parameters are set and described in covo.params


## input catalogues
Code searches for pairs of objects belonging to two different input catalogues.  

Input files should be ASCII and contain information on the position and  
components of two vectors for each object, i.e.  
x, y, z, vec_a_x, vec_a_y, vec_a_z, vec_b_x, vec_b_y, vec_b_z  

Both input catalogues must have the same format.  
Columns and column delimiter are set in parameter file.

