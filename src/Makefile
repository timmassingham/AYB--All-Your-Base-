CC = gcc
#CFLAGS = -Wall -ffast-math -O4 -funroll-loops -DNDEBUG -fomit-frame-pointer -std=gnu99 -msse2 -Wno-unused-function -fstrict-aliasing 
CFLAGS = -std=c99 -Wall -g -std=gnu99 -fstack-protector-all -D_FORTIFY_SOURCE=2 -Wno-unused-function
#CFLAGS = -O3 -pg -std=gnu99 
LD = ld
LDFLAGS = -lgfortran -lm -lz -lbz2 -lblas -llapack
INCFLAGS = 
#DEFINES = -D_GNU_SOURCE __linux
DEFINES = -DHAS_REALLOCF
objects = cif.o xio.o ayb.o matrix.o nuc.o estimate_lambda.o utility.o call_bases.o process_intensities.o estimateMPC.o statistics.o options.o main.o coordinates.o dnnls.o

all: ayb

ayb: $(objects)
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) -o ../bin/$@ $(LDFLAGS) $(objects) 


.c.o:
	$(CC) $(DEFINES) $(CFLAGS) $(INCFLAGS) -o $@ -c $<

.f.o:
	gfortran -m64 -O3 -o $@ -c $<

clean:
	rm -f ../bin/* *~ *.o  *.obj  *.exe *.lib

archive: clean
	cd ../.. && tar zcf ayb.tgz --exclude=.git --exclude=archive --exclude=test AYBc
