CXX=g++
CFLAGS=-I${GAMMALIB}/include/gammalib
LDFLAGS=-L${GAMMALIB}/lib -lgamma
DEPS=
OBJ=createfits.cpp

createfits: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS) $(LDFLAGS)
