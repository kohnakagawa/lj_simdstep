TARGET=a.out pair.out intrin.out

all: $(TARGET)

a.out: force.cpp
	icpc -O3 -xHOST -std=c++11 $< -o $@

pair.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DPAIR $< -o $@

intrin.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DINTRIN $< -o $@

clean:
	rm -f $(TARGET)

test: pair.out intrin.out
	./pair.out > pair.txt
	./intrin.out > intrin.txt
	diff pair.txt intrin.txt
