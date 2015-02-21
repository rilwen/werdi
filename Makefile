lapack=-L/opt/lapack/pathscale_3.0/3.1.1 -lpathfortran /opt/lapack/pathscale_3.0/3.1.1/core/lib64/liblapack.a /opt/lapack/pathscale_3.0/3.1.1/core/lib64/libblas.a

CC=pathcc
LIB_CCFLAGS=-Wall -O2 -std=c99 -march=core -I/home/awer/include
#LIB_CCFLAGS=-Wall -lefence -g -std=c99 -march=core -I/home/awer/include
#LIB_CCFLAGS=-Wall -pg -O2 -std=c99 -march=core -I/home/awer/include
TEST_CCFLAGS=-Wall -O2 -std=c99 -march=core -I/home/awer/include
#TEST_CCFLAGS=-Wall -lefence -g -std=c99 -march=core -I/home/awer/include
#TEST_CCFLAGS=-Wall -pg -O2 -std=c99 -march=core -I/home/awer/include
FC=pathf95
LIB_FCFLAGS=-Wall -O2 -march=core
HEADERS=doublecomplex.h constants.h hamiltonian.h tightbinding.h plane.h spinwaves.h
PREFIX=/home/awer
LIBDIR=$(PREFIX)/lib
INCLUDEDIR=$(PREFIX)/include

libwerdi=libwerdi.a

all: $(libwerdi)

install: $(libwerdi) werdi.h
	install -d $(INCLUDEDIR)/werdi
	install -m 644 *.h $(INCLUDEDIR)/werdi
	rm $(INCLUDEDIR)/werdi/werdi.h
	install -m 644 werdi.h $(INCLUDEDIR)
	install -m 644 $(libwerdi) $(LIBDIR)

tests: test_utils test_kane test_sixband test_secant test_berry test_layers test_spinwaves test_jancu test_lattice test_integrands test_find_EF test_anisotropy test_heig

$(libwerdi): utils.o fermi.o secant.o clapack.o kane.o clapack.o hamiltonian.o brillouin.o tightbinding.o sixband.o xi.o doublecomplex.o berry.o parabolic.o kp.o plane.o spinwaves.o anisotropy.o lattice.o fourband.o heig.o heigf.o $(HEADERS)
	$(CC) -ar -Wr,-r -o $@ $?
#	ar r $@ $?

%.o: %.c $(HEADERS)
	$(CC) $(LIB_CCFLAGS) -c $<

%.o: %.f
	$(FC) $(LIB_FCFLAGS) -c $<
	
test_utils: $(libwerdi) test_utils.c testing.o testing.h
	$(CC) $(TEST_CCFLAGS) test_utils.c testing.o -lm $(libwerdi) $(lapack) -o $@

test_kane: $(libwerdi) test_kane.c testing.o testing.h
	$(CC) $(TEST_CCFLAGS) test_kane.c testing.o -lm $(libwerdi) $(lapack) -o $@

test_sixband: $(libwerdi) test_sixband.c testing.o testing.h
	$(CC) $(TEST_CCFLAGS) test_sixband.c testing.o -lm $(libwerdi) $(lapack) -o $@

test_secant: $(libwerdi) test_secant.c testing.o testing.h
	$(CC) $(TEST_CCFLAGS) test_secant.c testing.o -lm $(libwerdi) $(lapack) -o $@
	
test_berry: $(libwerdi) test_berry.c testing.o testing.h
	$(CC) $(TEST_CCFLAGS) -L$(LIBDIR) test_berry.c testing.o -lm $(libwerdi) $(lapack) -ltb2 -o $@

test_layers: $(libwerdi) test_layers.c testing.o testing.h
	$(CC) $(TEST_CCFLAGS) -L$(LIBDIR) test_layers.c testing.o -lm $(libwerdi) $(lapack) -ltb2 -o $@
	
test_spinwaves: $(libwerdi) test_spinwaves.c testing.o testing.h
	$(CC) $(TEST_CCFLAGS) -L$(LIBDIR) test_spinwaves.c testing.o -lm $(libwerdi) $(lapack) -ltb2 -o $@

test_jancu: $(libwerdi) test_jancu.c testing.o testing.h
	$(CC) $(TEST_CCFLAGS) -L$(LIBDIR) test_jancu.c testing.o -lm $(libwerdi) $(lapack) -ltb2 -o $@

test_lattice: $(libwerdi) test_lattice.c testing.o testing.h
	$(CC) $(TEST_CCFLAGS) -L$(LIBDIR) test_lattice.c testing.o -lm $(libwerdi) $(lapack) -ltb2 -o $@

test_integrands: $(libwerdi) test_integrands.c testing.o testing.h
	$(CC) $(TEST_CCFLAGS) -L$(LIBDIR) test_integrands.c testing.o -lm $(libwerdi) $(lapack) -ltb2 -lcuba -o $@
	
test_find_EF: $(libwerdi) test_find_EF.c testing.o testing.h
	$(CC) $(TEST_CCFLAGS) -L$(LIBDIR) test_find_EF.c testing.o -lm $(libwerdi) $(lapack) -ltb2 -lcuba -o $@

test_anisotropy: $(libwerdi) test_anisotropy.c testing.o testing.h
	$(CC) $(TEST_CCFLAGS) -L$(LIBDIR) test_anisotropy.c testing.o -lm $(libwerdi) $(lapack) -ltb2 -lcuba -o $@

test_heig: $(libwerdi) test_heig.c testing.o testing.h
	$(CC) $(TEST_CCFLAGS) -L$(LIBDIR) test_heig.c testing.o -lm $(libwerdi) $(lapack) -ltb2 -lcuba -o $@

testing.o: testing.c doublecomplex.h
	$(CC) $(TEST_CCFLAGS) testing.c -c

clean:
	rm -rf *.o $(libwerdi)
