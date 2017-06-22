# -- Begin  force_sorted_intrin_mat_transp()
	.text
# mark_begin;
       .align    16,0x90
	.globl force_sorted_intrin_mat_transp()
# --- force_sorted_intrin_mat_transp()
force_sorted_intrin_mat_transp():
..B9.1:                         # Preds ..B9.0
..L487:         # optimization report
                # LOOP WAS FUSED
                # LOOP WAS COMPLETELY UNROLLED BY 1
                # LOOP WAS VECTORIZED
                # CILK PLUS ARRAY NOTATION LOOP
                # VECTORIZATION SPEEDUP COEFFECIENT 3.000000
                # LOOP HAS ONE VECTOR ITERATION
                # VECTOR LENGTH 4
                # MAIN VECTOR TYPE: 64-bits floating point
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
..___tag_value__Z30force_sorted_intrin_mat_transpv.474:
..L475:
                                                        #780.38
..LN2184:
	.loc    1  780  is_stmt 1
        pushq     %r12                                          #780.38
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
..LN2185:
        pushq     %r13                                          #780.38
	.cfi_def_cfa_offset 24
	.cfi_offset 13, -24
..LN2186:
        pushq     %r14                                          #780.38
	.cfi_def_cfa_offset 32
	.cfi_offset 14, -32
..LN2187:
        pushq     %r15                                          #780.38
	.cfi_def_cfa_offset 40
	.cfi_offset 15, -40
..LN2188:
        pushq     %rbx                                          #780.38
	.cfi_def_cfa_offset 48
	.cfi_offset 3, -48
..LN2189:
	.loc    1  782  prologue_end  is_stmt 1
        vmovupd   .L_2il0floatpacket.41(%rip), %ymm2            #782.21
..LN2190:
	.loc    1  783  is_stmt 1
        vmovupd   .L_2il0floatpacket.42(%rip), %ymm3            #783.21
..LN2191:
	.loc    1  784  is_stmt 1
        vmovupd   .L_2il0floatpacket.43(%rip), %ymm4            #784.21
..LN2192:
	.loc    1  781  is_stmt 1
        vxorpd    %ymm1, %ymm1, %ymm1                           #781.22
..LN2193:
	.loc    1  786  is_stmt 1
        xorl      %r8d, %r8d                                    #786.3
..LN2194:
	.loc    1  785  is_stmt 1
        movslq    particle_number(%rip), %r9                    #785.18
..LN2195:
	.loc    1  786  is_stmt 1
        xorl      %edi, %edi                                    #786.3
..LN2196:
        testq     %r9, %r9                                      #786.23
..LN2197:
        jle       ..B9.15       # Prob 10%                      #786.23
..LN2198:
                                # LOE rbx rbp rdi r8 r9 r12 r13 r14 r15 ymm1 ymm2 ymm3 ymm4
..B9.2:                         # Preds ..B9.1
..LN2199:
	.loc    1  855  is_stmt 1
        vmovsd    .L_2il0floatpacket.37(%rip), %xmm13           #855.16
..LN2200:
	.loc    1  857  is_stmt 1
        vmovsd    .L_2il0floatpacket.38(%rip), %xmm12           #857.59
..LN2201:
        vmovsd    .L_2il0floatpacket.39(%rip), %xmm11           #857.21
..LN2202:
        vmovsd    .L_2il0floatpacket.40(%rip), %xmm0            #857.33
..LN2203:
                                # LOE rbx rbp rdi r8 r9 r13 r14 r15 xmm0 xmm11 xmm12 xmm13 ymm1 ymm2 ymm3 ymm4
..B9.3:                         # Preds ..B9.13 ..B9.2
..L489:         # optimization report
                # LOOP WAS COMPLETELY UNROLLED BY 1
                # LOOP WAS VECTORIZED
                # CILK PLUS ARRAY NOTATION LOOP
                # VECTORIZATION SPEEDUP COEFFECIENT 4.667969
                # LOOP HAS ONE VECTOR ITERATION
                # VECTOR LENGTH 4
                # MAIN VECTOR TYPE: 64-bits floating point
..L488:         # optimization report
                # LOOP WITH USER VECTOR INTRINSICS
                # %s was not vectorized: inner loop was already vectorized
..LN2204:
	.loc    1  787  is_stmt 1
        vmovupd   q(%rdi), %ymm10                               #787.47
..LN2205:
	.loc    1  788  is_stmt 1
        vmovupd   p(%rdi), %ymm9                                #788.41
..LN2206:
	.loc    1  791  is_stmt 1
        xorl      %r12d, %r12d                                  #791.5
..LN2207:
	.loc    1  789  is_stmt 1
        movl      number_of_partners(,%r8,4), %ecx              #789.20
..LN2208:
	.loc    1  791  is_stmt 1
        movl      %ecx, %edx                                    #791.31
..LN2209:
        sarl      $1, %edx                                      #791.31
..LN2210:
        shrl      $30, %edx                                     #791.31
..LN2211:
        addl      %ecx, %edx                                    #791.31
..LN2212:
        sarl      $2, %edx                                      #791.31
..LN2213:
	.loc    1  790  is_stmt 1
        movslq    pointer(,%r8,4), %rsi                         #790.20
..LN2214:
	.loc    1  791  is_stmt 1
        lea       (,%rdx,4), %eax                               #791.36
..LN2215:
        testl     %eax, %eax                                    #791.36
..LN2216:
        jle       ..B9.7        # Prob 10%                      #791.36
..LN2217:
                                # LOE rbx rbp rsi rdi r8 r9 r12 r13 r14 r15 eax edx ecx xmm0 xmm11 xmm12 xmm13 ymm1 ymm2 ymm3 ymm4 ymm9 ymm10
..B9.4:                         # Preds ..B9.3
..LN2218:
        lea       3(,%rdx,4), %r10d                             #791.36
..LN2219:
        sarl      $1, %r10d                                     #791.36
..LN2220:
        shrl      $30, %r10d                                    #791.36
..LN2221:
        lea       3(%r10,%rdx,4), %r11d                         #791.36
..LN2222:
        sarl      $2, %r11d                                     #791.36
..LN2223:
	.loc    1  792  is_stmt 1
        lea       (,%rsi,4), %r10                               #792.23
..LN2224:
	.loc    1  791  is_stmt 1
        movslq    %r11d, %r11                                   #791.5
..LN2225:
                                # LOE rbp rsi rdi r8 r9 r10 r11 r12 eax edx ecx xmm0 ymm1 ymm2 ymm3 ymm4 ymm9 ymm10
..B9.5:                         # Preds ..B9.5 ..B9.4
..L506:         # optimization report
                # LOOP STMTS WERE REORDERED
                # LOOP WAS COMPLETELY UNROLLED BY 1
                # LOOP WAS VECTORIZED
                # CILK PLUS ARRAY NOTATION LOOP
                # VECTORIZATION SPEEDUP COEFFECIENT 5.000000
                # LOOP HAS ONE VECTOR ITERATION
                # VECTOR LENGTH 4
                # MAIN VECTOR TYPE: 64-bits floating point
..L505:         # optimization report
                # LOOP STMTS WERE REORDERED
                # LOOP WAS COMPLETELY UNROLLED BY 1
                # LOOP WAS VECTORIZED
                # CILK PLUS ARRAY NOTATION LOOP
                # VECTORIZATION SPEEDUP COEFFECIENT 5.000000
                # LOOP HAS ONE VECTOR ITERATION
                # VECTOR LENGTH 4
                # MAIN VECTOR TYPE: 64-bits floating point
..L504:         # optimization report
                # LOOP STMTS WERE REORDERED
                # LOOP WAS COMPLETELY UNROLLED BY 1
                # LOOP WAS VECTORIZED
                # CILK PLUS ARRAY NOTATION LOOP
                # VECTORIZATION SPEEDUP COEFFECIENT 5.000000
                # LOOP HAS ONE VECTOR ITERATION
                # VECTOR LENGTH 4
                # MAIN VECTOR TYPE: 64-bits floating point
..L503:         # optimization report
                # LOOP STMTS WERE REORDERED
                # LOOP WAS COMPLETELY UNROLLED BY 1
                # LOOP WAS VECTORIZED
                # CILK PLUS ARRAY NOTATION LOOP
                # VECTORIZATION SPEEDUP COEFFECIENT 5.046875
                # LOOP HAS ONE VECTOR ITERATION
                # VECTOR LENGTH 4
                # MAIN VECTOR TYPE: 64-bits floating point
..L502:         # optimization report
                # LOOP WAS COMPLETELY UNROLLED BY 1
                # LOOP WAS VECTORIZED
                # CILK PLUS ARRAY NOTATION LOOP
                # VECTORIZATION SPEEDUP COEFFECIENT 4.667969
                # LOOP HAS ONE VECTOR ITERATION
                # VECTOR LENGTH 4
                # MAIN VECTOR TYPE: 64-bits floating point
..L501:         # optimization report
                # LOOP WAS COMPLETELY UNROLLED BY 1
                # LOOP WAS VECTORIZED
                # CILK PLUS ARRAY NOTATION LOOP
                # VECTORIZATION SPEEDUP COEFFECIENT 4.667969
                # LOOP HAS ONE VECTOR ITERATION
                # VECTOR LENGTH 4
                # MAIN VECTOR TYPE: 64-bits floating point
..L500:         # optimization report
                # LOOP WAS COMPLETELY UNROLLED BY 1
                # LOOP WAS VECTORIZED
                # CILK PLUS ARRAY NOTATION LOOP
                # VECTORIZATION SPEEDUP COEFFECIENT 4.667969
                # LOOP HAS ONE VECTOR ITERATION
                # VECTOR LENGTH 4
                # MAIN VECTOR TYPE: 64-bits floating point
..L499:         # optimization report
                # LOOP WAS COMPLETELY UNROLLED BY 1
                # LOOP WAS VECTORIZED
                # CILK PLUS ARRAY NOTATION LOOP
                # VECTORIZATION SPEEDUP COEFFECIENT 4.667969
                # LOOP HAS ONE VECTOR ITERATION
                # VECTOR LENGTH 4
                # MAIN VECTOR TYPE: 64-bits floating point
..L498:         # optimization report
                # LOOP STMTS WERE REORDERED
                # LOOP WAS COMPLETELY UNROLLED BY 1
                # LOOP WAS VECTORIZED
                # CILK PLUS ARRAY NOTATION LOOP
                # VECTORIZATION SPEEDUP COEFFECIENT 4.351562
                # LOOP HAS ONE VECTOR ITERATION
                # VECTOR LENGTH 4
                # MAIN VECTOR TYPE: 64-bits floating point
..L497:         # optimization report
                # LOOP WAS COMPLETELY UNROLLED BY 1
                # LOOP WAS VECTORIZED
                # CILK PLUS ARRAY NOTATION LOOP
                # VECTORIZATION SPEEDUP COEFFECIENT 2.933594
                # LOOP HAS ONE VECTOR ITERATION
                # VECTOR LENGTH 4
                # MAIN VECTOR TYPE: 64-bits floating point
..L496:         # optimization report
                # LOOP WAS COMPLETELY UNROLLED BY 1
                # LOOP WAS VECTORIZED
                # CILK PLUS ARRAY NOTATION LOOP
                # VECTORIZATION SPEEDUP COEFFECIENT 2.933594
                # LOOP HAS ONE VECTOR ITERATION
                # VECTOR LENGTH 4
                # MAIN VECTOR TYPE: 64-bits floating point
..L495:         # optimization report
                # LOOP WAS COMPLETELY UNROLLED BY 1
                # LOOP WAS VECTORIZED
                # CILK PLUS ARRAY NOTATION LOOP
                # VECTORIZATION SPEEDUP COEFFECIENT 2.933594
                # LOOP HAS ONE VECTOR ITERATION
                # VECTOR LENGTH 4
                # MAIN VECTOR TYPE: 64-bits floating point
..L494:         # optimization report
                # LOOP WAS COMPLETELY UNROLLED BY 1
                # LOOP WAS VECTORIZED
                # CILK PLUS ARRAY NOTATION LOOP
                # VECTORIZATION SPEEDUP COEFFECIENT 2.933594
                # LOOP HAS ONE VECTOR ITERATION
                # VECTOR LENGTH 4
                # MAIN VECTOR TYPE: 64-bits floating point
..L493:         # optimization report
                # LOOP WAS COMPLETELY UNROLLED BY 1
                # LOOP WAS VECTORIZED
                # CILK PLUS ARRAY NOTATION LOOP
                # VECTORIZATION SPEEDUP COEFFECIENT 2.933594
                # LOOP HAS ONE VECTOR ITERATION
                # VECTOR LENGTH 4
                # MAIN VECTOR TYPE: 64-bits floating point
..L492:         # optimization report
                # LOOP WAS COMPLETELY UNROLLED BY 1
                # LOOP WAS VECTORIZED
                # CILK PLUS ARRAY NOTATION LOOP
                # VECTORIZATION SPEEDUP COEFFECIENT 2.933594
                # LOOP HAS ONE VECTOR ITERATION
                # VECTOR LENGTH 4
                # MAIN VECTOR TYPE: 64-bits floating point
..L491:         # optimization report
                # LOOP STMTS WERE REORDERED
                # LOOP WAS COMPLETELY UNROLLED BY 1
                # LOOP WAS VECTORIZED
                # CILK PLUS ARRAY NOTATION LOOP
                # VECTORIZATION SPEEDUP COEFFECIENT 4.949219
                # LOOP HAS ONE VECTOR ITERATION
                # VECTOR LENGTH 4
                # MAIN VECTOR TYPE: 64-bits floating point
..L490:         # optimization report
                # LOOP WITH USER VECTOR INTRINSICS
                # %s was not vectorized: inner loop was already vectorized
..LN2226:
	.loc    1  792  is_stmt 1
        movslq    sorted_list(%r10), %rbx                       #792.23
..LN2227:
	.loc    1  796  is_stmt 1
        movslq    4+sorted_list(%r10), %r13                     #796.23
..LN2228:
	.loc    1  800  is_stmt 1
        movslq    8+sorted_list(%r10), %r14                     #800.23
..LN2229:
	.loc    1  804  is_stmt 1
        movslq    12+sorted_list(%r10), %r15                    #804.23
..LN2230:
	.loc    1  793  is_stmt 1
        shlq      $5, %rbx                                      #793.45
..LN2231:
	.loc    1  797  is_stmt 1
        shlq      $5, %r13                                      #797.45
..LN2232:
	.loc    1  801  is_stmt 1
        shlq      $5, %r14                                      #801.45
..LN2233:
	.loc    1  805  is_stmt 1
        shlq      $5, %r15                                      #805.45
..LN2234:
	.loc    1  793  is_stmt 1
        vmovupd   q(%rbx), %ymm14                               #793.45
..LN2235:
	.loc    1  797  is_stmt 1
        vmovupd   q(%r13), %ymm6                                #797.45
..LN2236:
	.loc    1  801  is_stmt 1
        vmovupd   q(%r14), %ymm7                                #801.45
..LN2237:
	.loc    1  805  is_stmt 1
        vmovupd   q(%r15), %ymm8                                #805.45
..LN2238:
	.loc    1  798  is_stmt 1
        vsubpd    %ymm10, %ymm6, %ymm5                          #798.29
..LN2239:
	.loc    1  794  is_stmt 1
        vsubpd    %ymm10, %ymm14, %ymm6                         #794.29
..LN2240:
	.loc    1  806  is_stmt 1
        vsubpd    %ymm10, %ymm8, %ymm8                          #806.29
..LN2241:
	.loc    1  802  is_stmt 1
        vsubpd    %ymm10, %ymm7, %ymm7                          #802.29
..LN2242:
	.loc    1  808  is_stmt 1
        vunpcklpd %ymm5, %ymm6, %ymm13                          #808.19
..LN2243:
	.loc    1  809  is_stmt 1
        vunpckhpd %ymm5, %ymm6, %ymm12                          #809.19
..LN2244:
	.loc    1  810  is_stmt 1
        vunpcklpd %ymm8, %ymm7, %ymm15                          #810.19
..LN2245:
	.loc    1  811  is_stmt 1
        vunpckhpd %ymm8, %ymm7, %ymm14                          #811.19
..LN2246:
	.loc    1  813  is_stmt 1
        vperm2f128 $32, %ymm15, %ymm13, %ymm11                  #813.18
..LN2247:
	.loc    1  814  is_stmt 1
        vperm2f128 $32, %ymm14, %ymm12, %ymm14                  #814.18
..LN2248:
	.loc    1  815  is_stmt 1
        vperm2f128 $49, %ymm15, %ymm13, %ymm12                  #815.18
..LN2249:
	.loc    1  817  is_stmt 1
        vmulpd    %ymm14, %ymm14, %ymm15                        #817.36
..LN2250:
        vfmadd213pd %ymm15, %ymm11, %ymm11                      #817.36
..LN2251:
        vfmadd213pd %ymm11, %ymm12, %ymm12                      #817.48
..LN2252:
	.loc    1  818  is_stmt 1
        vmulpd    %ymm12, %ymm12, %ymm13                        #818.24
..LN2253:
        vmulpd    %ymm13, %ymm12, %ymm11                        #818.30
..LN2254:
	.loc    1  819  is_stmt 1
        vmulpd    %ymm11, %ymm11, %ymm15                        #819.47
..LN2255:
        vmovdqa   %ymm3, %ymm13                                 #819.32
..LN2256:
        vfmsub213pd %ymm4, %ymm11, %ymm13                       #819.32
..LN2257:
        vmulpd    %ymm15, %ymm12, %ymm11                        #819.53
..LN2258:
	.loc    1  820  is_stmt 1
        vsubpd    %ymm12, %ymm2, %ymm12                         #820.26
..LN2259:
	.loc    1  819  is_stmt 1
        vdivpd    %ymm11, %ymm13, %ymm14                        #819.53
..LN2260:
	.loc    1  821  is_stmt 1
        vblendvpd %ymm12, %ymm1, %ymm14, %ymm12                 #821.13
..LN2261:
	.loc    1  823  is_stmt 1
        vpermpd   $0, %ymm12, %ymm13                            #823.20
..LN2262:
	.loc    1  824  is_stmt 1
        vpermpd   $85, %ymm12, %ymm14                           #824.20
..LN2263:
	.loc    1  825  is_stmt 1
        vpermpd   $170, %ymm12, %ymm15                          #825.20
..LN2264:
	.loc    1  828  is_stmt 1
        vmovupd   p(%rbx), %ymm11                               #828.45
..LN2265:
	.loc    1  826  is_stmt 1
        vpermpd   $255, %ymm12, %ymm12                          #826.20
..LN2266:
	.loc    1  830  is_stmt 1
        vfnmadd231pd %ymm6, %ymm13, %ymm11                      #830.7
..LN2267:
	.loc    1  831  is_stmt 1
        vmovupd   %ymm11, p(%rbx)                               #831.33
..LN2268:
	.loc    1  833  is_stmt 1
        vmovupd   p(%r13), %ymm11                               #833.45
..LN2269:
	.loc    1  835  is_stmt 1
        vfnmadd231pd %ymm5, %ymm14, %ymm11                      #835.7
..LN2270:
	.loc    1  834  is_stmt 1
        vmulpd    %ymm14, %ymm5, %ymm5                          #834.7
..LN2271:
	.loc    1  836  is_stmt 1
        vmovupd   %ymm11, p(%r13)                               #836.33
..LN2272:
	.loc    1  834  is_stmt 1
        vfmadd213pd %ymm5, %ymm13, %ymm6                        #834.7
..LN2273:
	.loc    1  838  is_stmt 1
        vmovupd   p(%r14), %ymm5                                #838.45
..LN2274:
	.loc    1  834  is_stmt 1
        vaddpd    %ymm9, %ymm6, %ymm9                           #834.7
..LN2275:
	.loc    1  839  is_stmt 1
        vfmadd231pd %ymm7, %ymm15, %ymm9                        #839.7
..LN2276:
	.loc    1  840  is_stmt 1
        vfnmadd213pd %ymm5, %ymm15, %ymm7                       #840.7
..LN2277:
	.loc    1  841  is_stmt 1
        vmovupd   %ymm7, p(%r14)                                #841.33
..LN2278:
	.loc    1  843  is_stmt 1
        vmovupd   p(%r15), %ymm6                                #843.45
..LN2279:
	.loc    1  791  is_stmt 1
        incq      %r12                                          #791.5
..LN2280:
	.loc    1  844  is_stmt 1
        vfmadd231pd %ymm8, %ymm12, %ymm9                        #844.7
..LN2281:
	.loc    1  791  is_stmt 1
        addq      $16, %r10                                     #791.5
..LN2282:
	.loc    1  845  is_stmt 1
        vfnmadd213pd %ymm6, %ymm12, %ymm8                       #845.7
..LN2283:
	.loc    1  846  is_stmt 1
        vmovupd   %ymm8, p(%r15)                                #846.33
..LN2284:
	.loc    1  791  is_stmt 1
        cmpq      %r11, %r12                                    #791.5
..LN2285:
        jb        ..B9.5        # Prob 99%                      #791.5
..LN2286:
                                # LOE rbp rsi rdi r8 r9 r10 r11 r12 eax edx ecx xmm0 ymm1 ymm2 ymm3 ymm4 ymm9 ymm10
..B9.6:                         # Preds ..B9.5
..LN2287:
        vmovsd    .L_2il0floatpacket.39(%rip), %xmm11           #
..LN2288:
        vmovsd    .L_2il0floatpacket.38(%rip), %xmm12           #
..LN2289:
        vmovsd    .L_2il0floatpacket.37(%rip), %xmm13           #
..LN2290:
                                # LOE rbx rbp rsi rdi r8 r9 r13 r14 r15 eax edx ecx xmm0 xmm11 xmm12 xmm13 ymm1 ymm2 ymm3 ymm4 ymm9
..B9.7:                         # Preds ..B9.6 ..B9.3
..L507:         # optimization report
                # LOOP WAS COMPLETELY UNROLLED BY 1
                # LOOP WAS VECTORIZED
                # CILK PLUS ARRAY NOTATION LOOP
                # VECTORIZATION SPEEDUP COEFFECIENT 5.335938
                # LOOP HAS ONE VECTOR ITERATION
                # VECTOR LENGTH 4
                # MAIN VECTOR TYPE: 64-bits floating point
..LN2291:
	.loc    1  849  is_stmt 1
        xorl      %r10d, %r10d                                  #849.5
..LN2292:
        cmpl      %ecx, %eax                                    #849.36
..LN2293:
	.loc    1  848  is_stmt 1
        vmovupd   %ymm9, p(%rdi)                                #848.31
..LN2294:
	.loc    1  849  is_stmt 1
        jge       ..B9.13       # Prob 10%                      #849.36
..LN2295:
                                # LOE rbx rbp rsi rdi r8 r9 r10 r13 r14 r15 edx ecx xmm0 xmm11 xmm12 xmm13 ymm1 ymm2 ymm3 ymm4
..B9.8:                         # Preds ..B9.7
..LN2296:
	.loc    1  850  is_stmt 1
        movslq    %edx, %rdx                                    #850.21
..LN2297:
        movq      %rdx, %rax                                    #850.21
..LN2298:
        shlq      $4, %rax                                      #850.21
..LN2299:
	.loc    1  789  is_stmt 1
        shlq      $2, %rdx                                      #789.20
..LN2300:
	.loc    1  849  is_stmt 1
        movslq    %ecx, %rcx                                    #849.5
..LN2301:
	.loc    1  789  is_stmt 1
        negq      %rdx                                          #789.20
..LN2302:
	.loc    1  853  is_stmt 1
        vmovsd    16+q(%rdi), %xmm10                            #853.29
..LN2303:
	.loc    1  789  is_stmt 1
        addq      %rcx, %rdx                                    #789.20
..LN2304:
	.loc    1  852  is_stmt 1
        vmovsd    8+q(%rdi), %xmm9                              #852.29
..LN2305:
	.loc    1  850  is_stmt 1
        lea       (%rax,%rsi,4), %rax                           #850.21
..LN2306:
	.loc    1  851  is_stmt 1
        vmovsd    q(%rdi), %xmm5                                #851.29
..LN2307:
                                # LOE rax rdx rbx rbp rdi r8 r9 r10 r13 r14 r15 xmm0 xmm5 xmm9 xmm10 xmm11 xmm12 xmm13 ymm1 ymm2 ymm3 ymm4
..B9.9:                         # Preds ..B9.11 ..B9.8
..L508:         # optimization report
                # PEELED LOOP
                # LOOP STMTS WERE REORDERED
                # %s was not vectorized: vector dependence prevents vectorization%s
                # VECTOR TRIP COUNT IS ESTIMATED CONSTANT
..LN2308:
	.loc    1  850  is_stmt 1
        movslq    sorted_list(%rax,%r10,4), %rcx                #850.21
..LN2309:
	.loc    1  853  is_stmt 1
        shlq      $5, %rcx                                      #853.19
..LN2310:
	.loc    1  852  is_stmt 1
        vmovsd    8+q(%rcx), %xmm7                              #852.19
..LN2311:
	.loc    1  851  is_stmt 1
        vmovsd    q(%rcx), %xmm6                                #851.19
..LN2312:
	.loc    1  852  is_stmt 1
        vsubsd    %xmm9, %xmm7, %xmm7                           #852.29
..LN2313:
	.loc    1  851  is_stmt 1
        vsubsd    %xmm5, %xmm6, %xmm6                           #851.29
..LN2314:
	.loc    1  854  is_stmt 1
        vmulsd    %xmm7, %xmm7, %xmm14                          #854.35
..LN2315:
	.loc    1  853  is_stmt 1
        vmovsd    16+q(%rcx), %xmm8                             #853.19
..LN2316:
	.loc    1  854  is_stmt 1
        vfmadd231sd %xmm6, %xmm6, %xmm14                        #854.45
..LN2317:
	.loc    1  853  is_stmt 1
        vsubsd    %xmm10, %xmm8, %xmm8                          #853.29
..LN2318:
	.loc    1  854  is_stmt 1
        vfmadd231sd %xmm8, %xmm8, %xmm14                        #854.45
..LN2319:
	.loc    1  855  is_stmt 1
        vcomisd   %xmm13, %xmm14                                #855.16
..LN2320:
        ja        ..B9.11       # Prob 50%                      #855.16
..LN2321:
                                # LOE rax rdx rcx rbx rbp rdi r8 r9 r10 r13 r14 r15 xmm0 xmm5 xmm6 xmm7 xmm8 xmm9 xmm10 xmm11 xmm12 xmm13 xmm14 ymm1 ymm2 ymm3 ymm4
..B9.10:                        # Preds ..B9.9
..LN2322:
	.loc    1  856  is_stmt 1
        vmulsd    %xmm14, %xmm14, %xmm15                        #856.24
..LN2323:
        vmulsd    %xmm14, %xmm15, %xmm15                        #856.29
..LN2324:
	.loc    1  857  is_stmt 1
        vmulsd    %xmm15, %xmm14, %xmm14                        #857.47
..LN2325:
        vmulsd    %xmm15, %xmm14, %xmm14                        #857.52
..LN2326:
        vfmsub213sd %xmm0, %xmm11, %xmm15                       #857.59
..LN2327:
        vdivsd    %xmm14, %xmm15, %xmm15                        #857.52
..LN2328:
        vmulsd    %xmm15, %xmm12, %xmm14                        #857.59
..LN2329:
	.loc    1  858  is_stmt 1
        vmovsd    p(%rdi), %xmm15                               #858.7
..LN2330:
        vfmadd231sd %xmm6, %xmm14, %xmm15                       #858.7
..LN2331:
        vmovsd    %xmm15, p(%rdi)                               #858.7
..LN2332:
	.loc    1  859  is_stmt 1
        vmovsd    8+p(%rdi), %xmm15                             #859.7
..LN2333:
        vfmadd231sd %xmm14, %xmm7, %xmm15                       #859.7
..LN2334:
        vmovsd    %xmm15, 8+p(%rdi)                             #859.7
..LN2335:
	.loc    1  860  is_stmt 1
        vmovsd    16+p(%rdi), %xmm15                            #860.7
..LN2336:
        vfmadd231sd %xmm14, %xmm8, %xmm15                       #860.7
..LN2337:
        vmovsd    %xmm15, 16+p(%rdi)                            #860.7
..LN2338:
	.loc    1  861  is_stmt 1
        vmovsd    p(%rcx), %xmm15                               #861.7
..LN2339:
        vfnmadd213sd %xmm15, %xmm14, %xmm6                      #861.7
..LN2340:
        vmovsd    %xmm6, p(%rcx)                                #861.7
..LN2341:
	.loc    1  862  is_stmt 1
        vmovsd    8+p(%rcx), %xmm6                              #862.7
..LN2342:
        vfnmadd231sd %xmm14, %xmm7, %xmm6                       #862.7
..LN2343:
	.loc    1  863  is_stmt 1
        vmovsd    16+p(%rcx), %xmm7                             #863.7
..LN2344:
        vfnmadd213sd %xmm7, %xmm8, %xmm14                       #863.7
..LN2345:
	.loc    1  862  is_stmt 1
        vmovsd    %xmm6, 8+p(%rcx)                              #862.7
..LN2346:
	.loc    1  863  is_stmt 1
        vmovsd    %xmm14, 16+p(%rcx)                            #863.7
..LN2347:
                                # LOE rax rdx rbx rbp rdi r8 r9 r10 r13 r14 r15 xmm0 xmm5 xmm9 xmm10 xmm11 xmm12 xmm13 ymm1 ymm2 ymm3 ymm4
..B9.11:                        # Preds ..B9.9 ..B9.10
..LN2348:
	.loc    1  849  is_stmt 1
        incq      %r10                                          #849.5
..LN2349:
        cmpq      %rdx, %r10                                    #849.5
..LN2350:
        jb        ..B9.9        # Prob 99%                      #849.5
..LN2351:
                                # LOE rax rdx rbx rbp rdi r8 r9 r10 r13 r14 r15 xmm0 xmm5 xmm9 xmm10 xmm11 xmm12 xmm13 ymm1 ymm2 ymm3 ymm4
..B9.13:                        # Preds ..B9.11 ..B9.7
..LN2352:
	.loc    1  786  is_stmt 1
        incq      %r8                                           #786.3
..LN2353:
        addq      $32, %rdi                                     #786.3
..LN2354:
        cmpq      %r9, %r8                                      #786.3
..LN2355:
        jb        ..B9.3        # Prob 2%                       #786.3
..LN2356:
                                # LOE rbx rbp rdi r8 r9 r13 r14 r15 xmm0 xmm11 xmm12 xmm13 ymm1 ymm2 ymm3 ymm4
..B9.15:                        # Preds ..B9.13 ..B9.1
..LN2357:
	.loc    1  866  is_stmt 1
        vzeroupper                                              #866.1
	.cfi_restore 3
..LN2358:
	.loc    1  866  epilogue_begin  is_stmt 1
        popq      %rbx                                          #866.1
	.cfi_def_cfa_offset 40
	.cfi_restore 15
..LN2359:
        popq      %r15                                          #866.1
	.cfi_def_cfa_offset 32
	.cfi_restore 14
..LN2360:
        popq      %r14                                          #866.1
	.cfi_def_cfa_offset 24
	.cfi_restore 13
..LN2361:
        popq      %r13                                          #866.1
	.cfi_def_cfa_offset 16
	.cfi_restore 12
..LN2362:
        popq      %r12                                          #866.1
	.cfi_def_cfa_offset 8
..LN2363:
        ret                                                     #866.1
        .align    16,0x90
..LN2364:
                                # LOE
..LN2365:
	.cfi_endproc
# mark_end;
	.type	force_sorted_intrin_mat_transp(),@function
	.size	force_sorted_intrin_mat_transp(),.-force_sorted_intrin_mat_transp()
..LN_Z30force_sorted_intrin_mat_transpv.2366:
.LN_Z30force_sorted_intrin_mat_transpv:
	.data
# -- End  force_sorted_intrin_mat_transp()
