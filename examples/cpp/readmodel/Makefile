PROGRAM=readmodel
OBJ=readmodel.cpp
CFLAGS=`pkg-config --cflags gammalib`
LDFLAGS=`pkg-config --libs gammalib`

$(PROGRAM): $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS) $(LDFLAGS)