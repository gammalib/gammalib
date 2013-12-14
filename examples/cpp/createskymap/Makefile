CXX=g++
CFLAGS=-I${GAMMALIB}/include/gammalib
LDFLAGS=-L${GAMMALIB}/lib -lgamma
DEPS=
OBJ=createskymap.cpp

createskymap: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS) $(LDFLAGS)
