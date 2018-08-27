/* 
 * in0: in行列0 
 * in1, in行列1
 * out: 計算結果 
 * n: 辺の長さ 
 * pitch: パディング考慮した辺の長さ */

#include <util.h>
#include <arm_neon.h>
/* TYPE must be uint32_t */
void matmul_ver3(TYPE* in0, TYPE* in1, TYPE* out, size_t n, size_t pitch)
{
    /* i: 縦 j: 横 */
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            uint32x4_t l = vld1q_dup_u32(&in0[i * pitch + k]);
            for (int j = 0; j < n; j+=4) {
                //out[i * pitch + j] += in0[i * pitch + k] * in1[k * pitch + j];
                uint32x4_t o = vld1q_u32(&out[i * pitch + j]);
                uint32x4_t r = vld1q_u32(&in1[k * pitch + j]);
                o = vmlaq_u32(o, l, r);
                vst1q_u32(&out[i * pitch + j], o);
            }
        }
    }
}
