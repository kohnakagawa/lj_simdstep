#!/usr/bin/env python3

def is_comment(line):
    l = line.strip()
    return l[0] == '#'

def is_label(line):
    l = line.strip()
    return l.find(':') != -1

def is_directive(line):
    l = line.strip()
    return l[0] == '.'

def is_asminst(line):
    return not(is_comment(line) or is_label(line) or is_directive(line))

def is_lineinfo(line):
    return line.find(".loc") != -1

def get_src_lineidx(line):
    l = line.strip().split()
    return int(l[2])

src_f = "./force.cpp"

## for no gather/scatter 4
asm_f = "nogs_4_snip.asm"
code_region = (792, 846) # innermost loop
# code_region = (2041, 2169) # total
out_f = "nogs_4_insts.org"

src_line_num = range(code_region[0], code_region[1] + 1)
src_region_num = code_region[1] - code_region[0] + 1
insts_per_line = src_region_num * [0]

def is_line_in_region(src_line_idx):
    return src_line_idx >= code_region[0] and src_line_idx <= code_region[1]

with open(asm_f, "r") as f:
    asm_lines = f.readlines()
    num_asm_lines = len(asm_lines)
    cnt = 0
    src_line_idx = 0
    while cnt < num_asm_lines:
        l = asm_lines[cnt]
        cnt += 1
        if is_lineinfo(l):
            src_line_idx = get_src_lineidx(l)
        elif is_asminst(l):
            if is_line_in_region(src_line_idx):
                insts_per_line[src_line_idx - code_region[0]] += 1

with open(src_f, "r") as f_in, open(out_f, "w") as f_out:
    f_out.write("| # | src | # of insts | type |\n")
    for (i, l) in enumerate(f_in.readlines()):
        src_line_idx = i + 1
        if src_line_idx >= code_region[0] and src_line_idx <= code_region[1]:
            cnt_insts = insts_per_line[src_line_idx - code_region[0]]
            src_line = l.split("\n")[0]
            src_line = "{0:<100}".format(src_line)
            f_out.write("| %d | %s | %d | |\n" % (src_line_idx, src_line, cnt_insts))

print("# of total instructions")
print(sum(insts_per_line))
