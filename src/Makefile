LDFLAGS = -std=c++11 -O3
CC      = g++


OBJ = main.o parameters.o directory.o healpix.o catalogue.o  correlation.o
#------------------------------------------------------

covo: $(OBJ)
	$(CC) -o ../covo $(OBJ)


main.o: main.cpp
	${CC} -c main.cpp ${LDFLAGS}

parameters.o: parameters.cpp
	${CC} -c parameters.cpp ${LDFLAGS}
	
directory.o: directory.cpp
	${CC} -c directory.cpp ${LDFLAGS}

catalogue.o: catalogue.cpp
	${CC} -c catalogue.cpp ${LDFLAGS}

correlation.o: correlation.cpp
	${CC} -c correlation.cpp ${LDFLAGS}

healpix.o: healpix.cpp
	${CC} -c healpix.cpp ${LDFLAGS}
