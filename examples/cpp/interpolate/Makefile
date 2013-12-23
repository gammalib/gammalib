CXX=g++
CFLAGS=-I${GAMMALIB}/include/gammalib
LDFLAGS=-L${GAMMALIB}/lib -lgamma
DEPS=
OBJ=interpolate.cpp

interpolate: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS) $(LDFLAGS)
