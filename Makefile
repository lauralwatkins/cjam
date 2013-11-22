SHELL = /bin/sh

# ------------------------------------ #

# GSL libraries [path: -L/usr/local/lib is default]
GSLLIBS = -lgsl -fexternal-blas -framework Accelerate

# combine all libraries
LIBS = -lm $(GSLLIBS)

# compiler: -g produces debugging information, -Wall turns on all warnings
CC = gcc -g -Wall -ffast-math -O3 -fomit-frame-pointer

# compile options for compiling .c files to .o files
%/%.o: %/%.c
	$(CC) -fPIC -c $<

# ------------------------------------ #

INTERP = interp2dpol.o
INTERP := $(INTERP:%=interp/%)

JAM = jam_axi_rms_mgeint.o jam_axi_rms_mmt.o jam_axi_rms_wmmt.o \
	jam_axi_vel_losint.o jam_axi_vel_mgeint.o jam_axi_vel_mmt.o \
	jam_axi_vel_wmmt.o
JAM := $(JAM:%=jam/%)

MGE = mge_addbh.o mge_dens.o mge_deproject.o mge_qmed.o mge_read.o mge_surf.o
MGE := $(MGE:%=mge/%)

TOOLS = maximum.o median.o minimum.o range.o readcol.o sort_dbl.o where.o
TOOLS := $(TOOLS:%=tools/%)


cjam: $(INTERP) $(JAM) $(MGE) $(TOOLS) cjam.o cjam_main.o
	$(CC) $(INTERP) $(JAM) $(MGE) $(TOOLS) cjam.o cjam_main.o -o cjam $(LIBS) -L. -lpthread

clean: 
	rm *.o */*.o
