#ifndef __MATMUL_UTIL_H__
#define __MATMUL_UTIL_H__

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

#define TYPE uint32_t
typedef void (*CALLBACK_T)(TYPE*, TYPE*, TYPE*, size_t, size_t);

typedef struct {
    CALLBACK_T func;
    TYPE* in0;
    TYPE* in1;
    TYPE* out;
    size_t n;
    size_t pitch;
    double exec_time_sec;
    double gops;
} matmul_bench_t;

extern int verify_mat(matmul_bench_t* b, TYPE* ref);
extern void fill_rand_val(TYPE* val, size_t n, size_t pitch);
extern void fill_row_inc_val(TYPE* val, size_t n, size_t pitch);
extern void fill_col_inc_val(TYPE* val, size_t n, size_t pitch);
extern void fill_1(TYPE* val, size_t n, size_t pitch);
extern int bench_matmul(matmul_bench_t*);
extern void matmul_naive(TYPE* in0, TYPE* in1, TYPE* out, size_t n, size_t pitch);
extern void matmul_ver1(TYPE* in0, TYPE* in1, TYPE* out, size_t n, size_t pitch);
extern void matmul_ver2(TYPE* in0, TYPE* in1, TYPE* out, size_t n, size_t pitch);
extern void matmul_ver3(TYPE* in0, TYPE* in1, TYPE* out, size_t n, size_t pitch);
extern void matmul_ver4(TYPE* in0, TYPE* in1, TYPE* out, size_t n, size_t pitch);
#endif /* __MATMUL_UTIL_H__ */
