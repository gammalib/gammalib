CXX=g++
CFLAGS=-I${GAMMALIB}/include/gammalib
LDFLAGS=-L${GAMMALIB}/lib -lgamma
DEPS=
OBJ=readmodel.cpp

readmodel: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS) $(LDFLAGS)
