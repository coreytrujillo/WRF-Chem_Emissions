.SUFFIXES: .o .f90

FC     = ifort
F90    = $(FC)
LIBS   = -L$(NETCDF)/lib -lnetcdff -lnetcdf
INCLUDE_MODULES = -I$(NETCDF)/include
FFLAGS = -g

EXEC = fire_emis.exe 

OBJS = fire_utils.o\
       wrf_utils.o\
       glb_utils.o\
       fire_srf_types.o \
       fire_file.o \
       fire_emis.o

##dependencies
#$(OBJECTS) : makefile
%.o:             %.mod

${EXEC} :       ${OBJS}
		${F90} -o $@ ${FFLAGS} ${OBJS} ${LIBS} 

.f90.o:
		${F90} ${FFLAGS} -c ${INCLUDE_MODULES} $<

cleanup:
		rm -f ${OBJS} *.mod

clean:
		rm -f core ${EXEC} ${OBJS} *.mod
