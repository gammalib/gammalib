CXX=g++
CFLAGS=-I/usr/local/gamma/include/gammalib
LDFLAGS=-L/usr/local/gamma/lib -lgamma
DEPS=
OBJ=createskymap.cpp

createskymap: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS) $(LDFLAGS)
