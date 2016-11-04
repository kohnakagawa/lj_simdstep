TARGET=a.out pair.out s_intrin.out s_swp.out

all: $(TARGET)

a.out: force.cpp
	icpc -O3 -xHOST -std=c++11 $< -o $@

pair.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DPAIR $< -o $@

s_intrin.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DS_INTRIN $< -o $@

s_swp.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DS_SWP $< -o $@

clean:
	rm -f $(TARGET)

test: pair.out s_intrin.out
	./pair.out > pair.txt
	./s_intrin.out > s_intrin.txt
	diff pair.txt s_intrin.txt
