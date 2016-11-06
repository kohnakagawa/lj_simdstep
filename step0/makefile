PAIR=pair.out
SORTED=sorted.out s_swp.out s_intrin.out

TARGET=a.out $(PAIR) $(SORTED)

all: $(TARGET)

a.out: force.cpp
	icpc -O3 -xHOST -std=c++11 $< -o $@

pair.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DPAIR $< -o $@

sorted.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DSORTED $< -o $@

s_swp.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DS_SWP $< -o $@

s_intrin.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DS_INTRIN $< -o $@

clean:
	rm -f $(TARGET)

test: pair.out s_intrin.out
	./pair.out > pair.txt
	./s_intrin.out > s_intrin.txt
	diff pair.txt s_intrin.txt
