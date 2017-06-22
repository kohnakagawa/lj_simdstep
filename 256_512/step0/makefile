PAIR=pair.out
SORTED=sorted.out s_swp.out s_intrin.out

TARGET=a.out $(PAIR) $(SORTED)

CC=icpc
CPPFLAGS=-O3 -xHOST -std=c++11 -w3 -diag-disable:remark

all: $(TARGET)

a.out: force.cpp
	$(CC) $(CPPFLAGS) $< -o $@

pair.out: force.cpp
	$(CC) $(CPPFLAGS) -DPAIR $< -o $@

sorted.out: force.cpp
	$(CC) $(CPPFLAGS) -DSORTED $< -o $@

s_swp.out: force.cpp
	$(CC) $(CPPFLAGS) -DS_SWP $< -o $@

s_intrin.out: force.cpp
	$(CC) $(CPPFLAGS) -DS_INTRIN $< -o $@

clean:
	rm -f $(TARGET)

test: pair.out s_intrin.out
	./pair.out > pair.txt
	./s_intrin.out > s_intrin.txt
	diff pair.txt s_intrin.txt
