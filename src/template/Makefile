CXX=g++
CFLAGS=-I${GAMMALIB}/include/gammalib
LDFLAGS=-L${GAMMALIB}/lib -lgamma
DEPS=
OBJ=test_GClass.cpp GClass.cpp

test_GClass: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS) $(LDFLAGS)
