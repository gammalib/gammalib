CXX=g++
CFLAGS=-I${GAMMALIB}/include/gammalib
LDFLAGS=-L${GAMMALIB}/lib -lgamma
DEPS=
OBJ=createxml.cpp

createxml: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS) $(LDFLAGS)
