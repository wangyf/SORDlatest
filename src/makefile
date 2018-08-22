dest = /Users/yow004/Google_Drive_UCSD/code_under_development/SORDlatest/scripts/attenuation/run
# Compiler and Flags
#FC = gfortran
FC = mpif90

FLAGS = -O3

OBJS = surfnormals.o util.o frio.o bc.o globals.o diffcn.o diffnc.o \
       hourglass.o mpi.o arrays.o setup.o stats.o source.o checkpoint.o \
       fieldio.o parameters.o resample.o attenuation.o material.o timestep.o \
       inivolstress.o stress.o gridgen.o rupture.o acceleration.o sord.o

OBJSs = surfnormals.o util.o frio.o bc.o globals.o diffcn.o diffnc.o \
       hourglass.o serial.o arrays.o setup.o stats.o source.o checkpoint.o \
       fieldio.o parameters.o resample.o attenuation.o material.o timestep.o \
       inivolstress.o stress.o gridgen.o rupture.o acceleration.o sord.o

all:	sord-mO
	cp sord-mO $(dest)
   
sord-mO: $(OBJS)
	$(FC) $(FLAGS) $^ -o $@ 

sord-sO: $(OBJSs)
	$(FC) $(FLAGS) $^ -o $@ 
# independent
surfnormals.o:  surfnormals.f90
	$(FC) $(FLAGS) -c $^

util.o:         util.f90
	$(FC) $(FLAGS) -c $^

frio.o:         frio.f90
	$(FC) $(FLAGS) -c $^

bc.o:           bc.f90
	$(FC) $(FLAGS) -c $^

globals.o:      globals.f90
	$(FC) $(FLAGS) -c $^ 

diffcn.o:       diffcn.f90
	$(FC) $(FLAGS) -c $^

diffnc.o:       diffnc.f90
	$(FC) $(FLAGS) -c $^

hourglass.o:    hourglass.f90
	$(FC) $(FLAGS) -c $^

# dependent on m_frio
mpi.o:          mpi.f90
	$(FC) $(FLAGS) -c $^

serial.o:       serial.f90
	$(FC) $(FLAGS) -c $^
	
# dependent on globals.o
arrays.o:       arrays.f90
	$(FC) $(FLAGS) -c $^ 

# dependent on m_globals, m_collective (mpi.o), m_util
setup.o:        setup.f90
	$(FC) $(FLAGS) -c $^

stats.o:        stats.f90
	$(FC) $(FLAGS) -c $^

source.o:       source.f90
	$(FC) $(FLAGS) -c $^

# dependent on m_globals, m_collective (mpi.o), m_util, m_stats 
checkpoint.o:   checkpoint.f90
	$(FC) $(FLAGS) -c $^

# dependent on m_globals, m_collective (mpi.o), m_util, m_frio
fieldio.o:      fieldio.f90
	$(FC) $(FLAGS) -c $^

# dependent on m_globals, m_fieldio
parameters.o:   parameters.f90
	$(FC) $(FLAGS) -c $^

# dependent on m_globals, m_collective (mpi.o), m_util, m_bc
resample.o:     resample.f90
	$(FC) $(FLAGS) -c $^ 

# dependent on m_globals, m_collective (mpi.o), m_util, m_fieldio
material.o:     material.f90
	$(FC) $(FLAGS) -c $^

# dependent on m_globas, m_util, m_fieldio, m_stats 
timestep.o:     timestep.f90
	$(FC) $(FLAGS) -c $^

# dependent on m_globals, m_collective (mpi.o), m_util, m_bc, m_diffcn, m_fieldio
inivolstress.o: inivolstress.f90
	$(FC) $(FLAGS) -c $^ 

# dependent on m_globals, m_util, m_source, m_diffnc, m_fieldio, m_stats
stress.o:       stress.f90
	$(FC) $(FLAGS) -c $^

# dependent on m_globals, m_collective (mpi.o), m_util, m_bc, m_diffnc, m_fieldio, m_surfnormals
gridgen.o:      gridgen.f90
	$(FC) $(FLAGS) -c $^

# dependent on many above 
rupture.o:     rupture.f90
	$(FC) $(FLAGS) -c $^

acceleration.o: acceleration.f90
	$(FC) $(FLAGS) -c $^

attenuation.o:     attenuation.f90
	$(FC) $(FLAGS) -c $^

#diffcn.o:      diffcn.f90
#	$(FC) $(FLAGS) -c $^

# depdentent on all above
sord.o:         sord.f90
	$(FC) $(FLAGS) -c $^


# PHONY TARGETS 
clean:
	rm -rf $(OBJS) *.mod *.lst

cleanall:
	rm -rf *.o *.mod *.lst sord-* 

#cleanmore:
#	rm fort.* core.* CHK* SS* V* SRCT*
