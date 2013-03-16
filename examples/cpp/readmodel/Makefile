CXX=g++
CFLAGS=-I/usr/local/gamma/include/gammalib
LDFLAGS=-L/usr/local/gamma/lib -lgamma
DEPS=
OBJ=readmodel.cpp

%o: *.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

readmodel: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS) $(LDFLAGS)

