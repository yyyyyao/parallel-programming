/* 
 * in0: in行列0 
 * in1, in行列1
 * out: 計算結果 
 * n: 辺の長さ 
 * pitch: パディング考慮した辺の長さ */

#include <util.h>
#include <arm_neon.h>
/* TYPE must be uint32_t */
/* use register blocking */

#define NEON_OP(I, J) \
    vout ## I ## _ ## J = vmlaq_u32(vout ## I ## _ ## J, lik ## I ## _4, vr ## J);

#define NEON_K4(K) \
    { \
        const uint32_t* inL00 = inL0; \
         \
        lik0_4 = vld1q_dup_u32(inL00); \
        inL00 += pitch; \
        vr0 = vld1q_u32(&inR0[0]); \
         \
        NEON_OP(0, 0); \
        lik1_4 = vld1q_dup_u32(inL00); \
        inL00 += pitch; \
        NEON_OP(1, 0); \
        lik2_4 = vld1q_dup_u32(inL00); \
        NEON_OP(2, 0); \
         \
        vr1 = vld1q_u32(&inR0[4]); \
        NEON_OP(0, 1); \
        NEON_OP(1, 1); \
        NEON_OP(2, 1); \
         \
        vr2 = vld1q_u32(&inR0[8]); \
        NEON_OP(0, 2); \
        NEON_OP(1, 2); \
        NEON_OP(2, 2); \
         \
        vr3 = vld1q_u32(&inR0[12]); \
        NEON_OP(0, 3); \
        NEON_OP(1, 3); \
        NEON_OP(2, 3); \
        inL0++; \
        inR0 += pitch; \
    }

static void loop(int i, int j, int k, int bi,
        TYPE* out, TYPE* in0, TYPE* in1, size_t n, size_t pitch)
{
    int i0 = i + bi;
    TYPE* outp0 = &out[i0 * pitch + j];
    TYPE* outp1 = outp0 + (pitch * 1);
    TYPE* outp2 = outp0 + (pitch * 2);

    uint32x4_t* outp4_0 = (uint32x4_t*)outp0;
    uint32x4_t* outp4_1 = (uint32x4_t*)outp1;
    uint32x4_t* outp4_2 = (uint32x4_t*)outp2;

    uint32x4_t vout0_0 = vmovq_n_u32(0);
    uint32x4_t vout0_1 = vmovq_n_u32(0);
    uint32x4_t vout0_2 = vmovq_n_u32(0);
    uint32x4_t vout0_3 = vmovq_n_u32(0);

    uint32x4_t vout1_0 = vmovq_n_u32(0);
    uint32x4_t vout1_1 = vmovq_n_u32(0);
    uint32x4_t vout1_2 = vmovq_n_u32(0);
    uint32x4_t vout1_3 = vmovq_n_u32(0);

    uint32x4_t vout2_0 = vmovq_n_u32(0);
    uint32x4_t vout2_1 = vmovq_n_u32(0);
    uint32x4_t vout2_2 = vmovq_n_u32(0);
    uint32x4_t vout2_3 = vmovq_n_u32(0);

    /* calculation */
    uint32x4_t lik0_4, lik1_4, lik2_4;
    uint32x4_t vr0, vr1, vr2, vr3;

    const uint32_t *inL0 = &in0[i0 * pitch + k];
    const uint32_t *inR0 = &in1[k * pitch + j];

    for(int i = 0; i < 32; i++) {
#if 1
        NEON_K4(0);
        NEON_K4(1);
        NEON_K4(2);
        NEON_K4(3);
#else
        const uint32_t* inL00 = inL0;
        
        lik0_4 = vld1q_dup_u32(inL00);
        inL00 += pitch;
        vr0 = vld1q_u32(&inR0[0]);
        
        NEON_OP(0, 0);
        lik1_4 = vld1q_dup_u32(inL00);
        inL00 += pitch;
        NEON_OP(1, 0);
        lik2_4 = vld1q_dup_u32(inL00);
        NEON_OP(2, 0);
        
        vr1 = vld1q_u32(&inR0[4]);
        NEON_OP(0, 1);
        NEON_OP(1, 1);
        NEON_OP(2, 1);
        
        vr2 = vld1q_u32(&inR0[8]);
        NEON_OP(0, 2);
        NEON_OP(1, 2);
        NEON_OP(2, 2);
        
        vr3 = vld1q_u32(&inR0[12]);
        NEON_OP(0, 3);
        NEON_OP(1, 3);
        NEON_OP(2, 3);
        inL0++;
        inR0 += pitch;

        inL00 = inL0;

        lik0_4 = vld1q_dup_u32(inL00);
        inL00 += pitch;
        vr0 = vld1q_u32(&inR0[0]);
        
        NEON_OP(0, 0);
        lik1_4 = vld1q_dup_u32(inL00);
        inL00 += pitch;
        NEON_OP(1, 0);
        lik2_4 = vld1q_dup_u32(inL00);
        NEON_OP(2, 0);
        
        vr1 = vld1q_u32(&inR0[4]);
        NEON_OP(0, 1);
        NEON_OP(1, 1);
        NEON_OP(2, 1);
        
        vr2 = vld1q_u32(&inR0[8]);
        NEON_OP(0, 2);
        NEON_OP(1, 2);
        NEON_OP(2, 2);
        
        vr3 = vld1q_u32(&inR0[12]);
        NEON_OP(0, 3);
        NEON_OP(1, 3);
        NEON_OP(2, 3);
        inL0++;
        inR0 += pitch;
#endif
    }

    if (k == 0) {
        outp4_0[0] = vout0_0;
        outp4_0[1] = vout0_1;
        outp4_0[2] = vout0_2;
        outp4_0[3] = vout0_3;

        outp4_1[0] = vout1_0;
        outp4_1[1] = vout1_1;
        outp4_1[2] = vout1_2;
        outp4_1[3] = vout1_3;

        outp4_2[0] = vout2_0;
        outp4_2[1] = vout2_1;
        outp4_2[2] = vout2_2;
        outp4_2[3] = vout2_3;
    } else {
        outp4_0[0] = vaddq_u32(outp4_0[0], vout0_0);
        outp4_0[1] = vaddq_u32(outp4_0[1], vout0_1);
        outp4_0[2] = vaddq_u32(outp4_0[2], vout0_2);
        outp4_0[3] = vaddq_u32(outp4_0[3], vout0_3);

        outp4_1[0] = vaddq_u32(outp4_1[0], vout1_0);
        outp4_1[1] = vaddq_u32(outp4_1[1], vout1_1);
        outp4_1[2] = vaddq_u32(outp4_1[2], vout1_2);
        outp4_1[3] = vaddq_u32(outp4_1[3], vout1_3);

        outp4_2[0] = vaddq_u32(outp4_2[0], vout2_0);
        outp4_2[1] = vaddq_u32(outp4_2[1], vout2_1);
        outp4_2[2] = vaddq_u32(outp4_2[2], vout2_2);
        outp4_2[3] = vaddq_u32(outp4_2[3], vout2_3);
    }

}
void matmul_ver4(TYPE* in0, TYPE* in1, TYPE* out, size_t n, size_t pitch)
{
    const int block_size_i = 48;
    const int block_size_j = 32;
    const int block_size_k = 64;
    const int i_inc = 3;
    /* i: 縦 j: 横 */
    for (int i = 0; i < n; i+=block_size_i) {
        for (int j = 0; j < n; j+=block_size_j) {
            for (int k = 0; k < n; k+=block_size_k) {
                for(int bi = 0; bi < block_size_i; bi+= i_inc) {
                    loop(i, j, k, bi, out, in0, in1, n, pitch);
                }
            }
        }
    }
}
