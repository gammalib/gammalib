CXX=g++
CFLAGS=-I${GAMMALIB}/include/gammalib
LDFLAGS=-L${GAMMALIB}/lib -lgamma
DEPS=
OBJ=optimize.cpp

optimize: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS) $(LDFLAGS)
