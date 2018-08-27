/* 
 * in0: in行列0 
 * in1, in行列1
 * out: 計算結果 
 * n: 辺の長さ 
 * pitch: パディング考慮した辺の長さ */

#include <util.h>

void matmul_naive(TYPE* in0, TYPE* in1, TYPE* out, size_t n, size_t pitch)
{
    /* i: 縦 j: 横 */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            TYPE val = 0;
            for (int k = 0; k < n; k++) {
                val += in0[i * pitch + k] * in1[k * pitch + j];
            }
            out[i * pitch + j] = val;
        }
    }
}
