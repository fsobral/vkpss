all: lib

install:
	mkdir -p $(LIB)
	mv libsvr.a $(LIB)

lib: cinterface.o trdf.o TRDF.o svrmod.o
	ar crv libsvr.a $^

trdf.o: trdf.f90 svrmod.o
	$(FC) -c -J$(LIB) $< -o $@ 

cinterface.o: cinterface.f90
	$(FC) -c -J$(LIB) $< -o $@ 

svrmod.o: svrmod.f90
	$(FC) -c -J$(LIB) $< -o $@ 

TRDF.o: TRDF.f	
	$(FC) -c -J$(LIB) $^ -o $@

%.o: %.f90
	$(FC) -c -J$(LIB) $^ -o $@

clean:
	rm -vf *.o *.mod *~

.PHONY: lib clean all install
