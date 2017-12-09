AOS = force.out
ASM = force.s

TARGET= $(AOS) $(ASM)

OPT_FLAGS=-O3 -xHOST -no-prec-div -ipo

all: $(TARGET)

.SUFFIXES:
.SUFFIXES: .cpp .s
.cpp.s:
	icpc $(OPT_FLAGS) -std=c++11 -masm=intel -S $< -o $@

.SUFFIXES: .cpp .out
.cpp.out:
	icpc $(OPT_FLAGS) -std=c++11 $< -o $@

clean:
	rm -f $(TARGET) $(ASM)
