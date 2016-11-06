PAIR=pair.out p_swp.out
SORTED=sorted.out s_swp.out s_intrin.out

TARGET=a.out $(PAIR) $(SORTED)

all: $(TARGET)

a.out: force.cpp
	icpc -O3 -xHOST -std=c++11 $< -o $@

pair.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DPAIR $< -o $@

p_swp.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DP_SWP $< -o $@

sorted.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DSORTED $< -o $@

s_swp.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DS_SWP $< -o $@

s_intrin.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DS_INTRIN $< -o $@

clean:
	rm -f $(TARGET)

test: p_swp.out s_intrin.out
	./p_swp.out > p_swp.txt
	./s_intrin.out > s_intrin.txt
	diff p_swp.txt s_intrin.txt
