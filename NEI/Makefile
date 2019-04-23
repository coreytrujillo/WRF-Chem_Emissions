.SUFFIXES: .F

FC     = ifort
F90    = $(FC)
# LIBS   = -L$(NETCDF)/lib -lnetcdff -lnetcdf
# INCLUDE_MODULES = -I$(NETCDF)/include
FFLAGS = -axP -free -convert big_endian -o

BASE = emiss_v04_Rimfire_4km
EXEC = ${BASE}.exe 
OBJ = ${BASE}.F

${EXEC}:		${OBJ}
		${FC} ${FFLAGS} $@ $<

cleanup:
		rm -f ${OBJS} *.mod 

clean:
		rm -f core ${EXEC} ${OBJS} *.mod
