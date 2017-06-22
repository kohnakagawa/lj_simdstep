SORTED=sorted.out s_intrin_gs.out s_intrin_nogs.out
ASM=force.s
TARGET=a.out $(SORTED) $(ASM)

CC=icpc
CPPFLAGS=-O3 -xHOST -std=c++11 -w3 -diag-disable:remark

all: $(TARGET)

.SUFFIXES:
.SUFFIXES: .cpp .s
.cpp.s:
	$(CC) $(CPPFLAGS) -g -S $< -o $@

a.out: force.cpp
	$(CC) $(CPPFLAGS) $< -o $@

sorted.out: force.cpp
	$(CC) $(CPPFLAGS) -DSORTED $< -o $@

s_intrin_gs.out: force.cpp
	$(CC) $(CPPFLAGS) -DS_INTRIN_GS $< -o $@

s_intrin_nogs.out: force.cpp
	$(CC) $(CPPFLAGS) -DS_INTRIN_NOGS $< -o $@

clean:
	rm -f $(TARGET)
