TARGET=a.out pair.out s_intrin.out s_swp.out p_swp.out p_intrin.out

all: $(TARGET)

a.out: force.cpp
	icpc -O3 -xHOST -std=c++11 $< -o $@

pair.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DPAIR $< -o $@

p_swp.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DP_SWP $< -o $@

p_intrin.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DP_INTRIN $< -o $@

s_intrin.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DS_INTRIN $< -o $@

s_intrin_swp.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DS_SWP_INTRIN $< -o $@

s_swp.out: force.cpp
	icpc -O3 -xHOST -std=c++11 -DS_SWP $< -o $@

clean:
	rm -f $(TARGET)

test: s_intrin.out s_intrin_swp.out
	./s_intrin_swp.out > s_intrin_swp.txt
	./s_intrin.out > s_intrin.txt
	diff s_intrin_swp.txt s_intrin.txt
