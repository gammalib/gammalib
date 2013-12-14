CXX=g++
CFLAGS=-I${GAMMALIB}/include/gammalib
LDFLAGS=-L${GAMMALIB}/lib -lgamma
DEPS=
OBJ=numerics.cpp

numerics: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS) $(LDFLAGS)
