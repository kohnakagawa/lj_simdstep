PAIR=pair.out p_swp.out p_swp_intrin.out
SORTED=sorted.out s_swp.out s_intrin.out s_swp_intrin.out s_swp_intrin_m_trans.out

TARGET=a.out $(PAIR) $(SORTED) knl.out

CC=icpc
CPPFLAGS=-O3 -xHOST -std=c++11 -w3 -diag-disable:remark

all: $(TARGET)

a.out: force.cpp
	$(CC) $(CPPFLAGS) $< -o $@

knl.out: force.cpp
	$(CC) $(CPPFLAGS) -DKNL $< -o $@

knl_swp.out: force.cpp
	$(CC) $(CPPFLAGS) -DKNL_SWP $< -o $@

knl_gs.out: force.cpp
	$(CC) $(CPPFLAGS) -DKNL_GS $< -o $@

pair.out: force.cpp
	$(CC) $(CPPFLAGS) -DPAIR $< -o $@

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

test: s_intrin.out s_swp_intrin_m_trans.out
	./s_swp_intrin_m_trans.out > s_swp_intrin_m_trans.txt
	./s_intrin.out > s_intrin.txt
	diff s_swp_intrin_m_trans.txt s_intrin.txt

testz: knl.out sorted.out
	./sorted.out > sorted.txt
	./knl.out > knl.txt
	diff sorted.txt knl.txt

-include makefile.opt
