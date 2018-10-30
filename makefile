
CC = g++
CFLAGS = -lm -llapack -I.
DEPS = action.h evolve_fields.h force.h kinetic_energy.h matrixmodel.h measure.h my_gen.h myrandom.h obs.h read_in.h read_param.h setup.h update.h utilities.h write_out.h
OBJ = matrixmodel.o action.o evolve_fields.o force.o kinetic_energy.o matrixmodel.o measure.o my_gen.o myrandom.o obs.o read_in.o read_param.o setup.o update.o utilities.o write_out.o

%.o : %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)							

#I think the ifndef in simulation.h will become effective when doing the $< part because $(DEPS) will be used even with those *.cpp which don't need $(DEPS) and also those *.cpp where its already included.
#Explore the $< part properly.

output : $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)


clean :
	rm *.o output 
	clear

