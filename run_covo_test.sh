# option 0:
# all parameters are specified by the parameter file
#./covo covo.params_test catalogues/mice_halos.csv

# option 1:
# catalogue_1 is input argument, catalogue_2 is copied from catalogue_1
./covo covo.params_test ./catalogues/random_1.csv

# option 2:
# catalogue_1 and catalogue_2 are input arguments
#./covo covo.params ./catalogues/random_1.csv ./catalogues/random_2.csv

# option 3:
# suffix for output file as last argument
#./covo covo.params ./catalogues/random_1.csv ./catalogues/random_2.csv test
