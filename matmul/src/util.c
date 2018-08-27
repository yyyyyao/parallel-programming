#include <util.h>

double
static matmul_bench_calc_gops(unsigned int n, double sec)
{
    /* 要素数(n * n) * 各要素の計算量((乗算 + 加算) * n) 
     * = n * n * 2 * n */
    return (n*2.0*n*n) / (sec*1000.0*1000.0*1000.0);
}

static double matmul_bench_sec(void)
{
    struct timespec ts;

#ifdef CLOCK_MONOTONIC_RAW
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
#else
    clock_gettime(CLOCK_MONOTONIC, &ts);
#endif

    return (ts.tv_sec) + (ts.tv_nsec / (1000.0*1000.0*1000.0));
}

static void exec_func(matmul_bench_t* b) {
    b->func(b->in0, b->in1, b->out, b->n, b->pitch);
}

int bench_matmul(matmul_bench_t* bench) {
    double start, end;

    /* 設定されたコンフィグでmatmulを実行. */
    start = matmul_bench_sec();

    exec_func(bench);
    end = matmul_bench_sec();

    /* 結果を構造体の中に格納する */
    bench->exec_time_sec = end - start;
    bench->gops =
        matmul_bench_calc_gops(bench->n, bench->exec_time_sec);
}

void fill_rand_val(TYPE* val, size_t n, size_t pitch)
{
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            val[i * pitch + j] = rand() % 255;
        }
    }
}

void fill_col_inc_val(TYPE* val, size_t n, size_t pitch)
{
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            val[i * pitch + j] = i * n + j;
        }
    }
}

void fill_1(TYPE* val, size_t n, size_t pitch)
{
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            val[i * pitch + j] = 1;
        }
    }
}

void fill_row_inc_val(TYPE* val, size_t n, size_t pitch)
{
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            val[i * pitch + j] = j * n + i;
        }
    }
}
int verify_mat(matmul_bench_t* b, TYPE* ref)
{
    size_t pitch = b->pitch;
    for(int i = 0; i < b->n; i++) {
        for(int j = 0; j < b->n; j++) {
            int idx = i * pitch + j;
            if(b->out[idx] != ref[idx]) {
                printf("diff: %d %d out:%d ref:%d\n", i, j, 
                    b->out[idx],
                    ref[idx]);
                return -1;
           /*  }else {printf("ok\n"); */
            }
        }
    }
    return 0;
}
