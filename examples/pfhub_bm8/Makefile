# If you have a different compiler then you can 
# specify this by override, e.g.: make CC='pgcc'
CC = g++
CCOPTS = -g
PGOPTS = 
LIBS = -lm

incdir = include

# If you have gnuplot pipes on your system, change this
# value to 'yes' and change the options in the source
# code to get graphical output.
GNUPLOT_I = no

# If different compilers have different options, then
# those can be handled here.
ifeq ($(CC),gcc)
	CCOPTIONS = $(CCOPTS)
endif

ifeq ($(CC),pgcc)
	CCOPTIONS = $(PGOPTS)
endif

ifeq ($(CC),g++)
	CCOPTIONS = 
endif

ifeq ($(GNUPLOT_I),yes)
	OBJECT1 = MMSPNucleation.o   
else
	OBJECT1 = MMSPNucleation.o
endif

#Program dependencies below here.
all :  MMSPNucleation
.PHONY : all

MMSPNucleation : $(OBJECT1)
	$(CC) $(CCOPTIONS) -o $@ MMSPNucleation.o $(LIBS)


# Object defs below here.

MMSPNucleation.o : MMSPNucleation.cpp
	$(CC) $(CCOPTIONS) -I$(incdir) -c $<

gnuplot_i.hpp.gch : gnuplot_i.hpp
	$(CC) $(CCOPTIONS) -c gnuplot_i.hpp

clean :
	rm -rf 1stDiffusion MMSPDiffusion MMSPDiffusion2D *.o *~ *.gch

