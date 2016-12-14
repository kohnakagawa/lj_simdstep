PAIR=pair.out p_swp.out p_swp_intrin.out
SORTED=sorted.out s_swp.out s_intrin.out s_swp_intrin.out s_swp_intrin_m_trans.out pair_wsort.out

TARGET=a.out $(PAIR) $(SORTED)

CC=icpc
CPPFLAGS=-O3 -xHOST -std=c++11 -w3 -diag-disable:remark

all: $(TARGET)

a.out: force.cpp
	$(CC) $(CPPFLAGS) $< -o $@

pair.out: force.cpp
	$(CC) $(CPPFLAGS) -DPAIR $< -o $@

pair_wsort.out: force.cpp
	$(CC) $(CPPFLAGS) -DPAIR -DWITH_SORT $< -o $@

p_swp.out: force.cpp
	$(CC) $(CPPFLAGS) -DP_SWP $< -o $@

p_swp_intrin.out: force.cpp
	$(CC) $(CPPFLAGS) -DP_SWP_INTRIN $< -o $@

sorted.out: force.cpp
	$(CC) $(CPPFLAGS) -DSORTED $< -o $@

s_swp.out: force.cpp
	$(CC) $(CPPFLAGS) -DS_SWP $< -o $@

s_intrin.out: force.cpp
	$(CC) $(CPPFLAGS) -DS_INTRIN $< -o $@

s_swp_intrin.out: force.cpp
	$(CC) $(CPPFLAGS) -DS_SWP_INTRIN $< -o $@

s_swp_intrin_m_trans.out: force.cpp
	$(CC) $(CPPFLAGS) -DS_SWP_INTRIN_M_TRANS $< -o $@

clean:
	rm -f $(TARGET)

clear:
	rm -f .pair.dat .pair_sorted.dat

test: s_intrin.out s_swp_intrin.out
	./s_swp_intrin.out > s_swp_intrin.txt
	./s_intrin.out > s_intrin.txt
	diff s_swp_intrin.txt s_intrin.txt
