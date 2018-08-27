#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"

#define VERIFY_RESULT

int main(int argc, char* argv[])
{
    matmul_bench_t bench;
    TYPE* ref;
    //size_t len_of_side[] = {256, 384, 512, 640, 768, 1024};
    //size_t len_of_side[] = {256, 384};
    size_t len_of_side[] = {128};
    CALLBACK_T funcs[] = {matmul_ver4};

    for(int i = 0; i < sizeof(len_of_side) / sizeof(size_t); i++) {
        size_t len = len_of_side[i];
        /* for cache */
        size_t alloc_bytes = sizeof(TYPE) * len * (len + 1);

        /* alloc */
        bench.in0 = (TYPE*)aligned_alloc(sizeof(TYPE), alloc_bytes);
        bench.in1 = (TYPE*)aligned_alloc(sizeof(TYPE), alloc_bytes);
        bench.out = (TYPE*)aligned_alloc(sizeof(TYPE), alloc_bytes);
        bench.n   = len;
        bench.pitch = len + 1;

#ifdef VERIFY_RESULT
        /* for ref */
        ref = (TYPE*)aligned_alloc(sizeof(TYPE), alloc_bytes);
#endif

        /* initialize in0, in1 */
        //fill_rand_val(bench.in0, bench.n, bench.pitch);
        //fill_rand_val(bench.in1, bench.n, bench.pitch);

        //fill_col_inc_val(bench.in0, bench.n, bench.pitch);
        //fill_row_inc_val(bench.in1, bench.n, bench.pitch);
        
        fill_1(bench.in0, bench.n, bench.pitch);
        fill_1(bench.in1, bench.n, bench.pitch);

        /* create ref */
#ifdef VERIFY_RESULT
        memset(bench.out, 0x0, alloc_bytes);
        bench.func = matmul_naive;
        bench_matmul(&bench);
        memcpy(ref, bench.out, alloc_bytes);
#endif

        for(int j = 0; j < sizeof(funcs) / sizeof(CALLBACK_T); j++) {
            /* initialize out */
            memset(bench.out, 0x0, alloc_bytes);

            bench.func = funcs[j];

            bench_matmul(&bench);
#ifdef VERIFY_RESULT
            if (verify_mat(&bench, ref) != 0) {
                printf("bench_matmul result is incorrect!\n!");
                return -1;
            }
#endif
            printf("%ldx%ld pitch:%ld %f[GOPS] %f[sec]\n",
                    bench.n, bench.n, bench.pitch, bench.gops, bench.exec_time_sec);
        }
        /* free */
        free(bench.in0);
        free(bench.in1);
        free(bench.out);
#ifdef VERIFY_RESULT
        free(ref);
#endif
    }

    return 0;
}
