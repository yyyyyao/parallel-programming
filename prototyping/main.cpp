#include <stdio.h>
#include <malloc.h>
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <tmmintrin.h>

typedef union
{
	float a[4];
	__m128 v;
}SSEFloat;

void printDat(float*, int, int);
void printRefLineVertical(float*, int, int);
void printXmm(__m128 data) {
	SSEFloat s;
	s.v = data;
	printf("0:%d 1:%d 2:%d 3:%d\n",
		(int)s.a[0], (int)s.a[1], (int)s.a[2], (int)s.a[3]);
}

int main(void) {

	int width = 7;
	int height = 7;
	float* dat;
	dat = (float*)malloc(sizeof(float) * width * height);

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			dat[i * width + j] = (float)(i * width + j);
		}
	}

	__m128 tmm0, tmm1, tmm2, tmm3, tmm4, tmm5, tmm6, tmm7;
	float aaa, bbb, ccc, ddd, eee, fff, ggg, hhh, iii, jjj, kkk, lll, mmm, nnn,
		ooo, ppp, qqq, rrr;
	SSEFloat u;
	float* buf;
	buf = (float*)_aligned_malloc(sizeof(float) * 64, 16);

	/* 1-4 */
	buf[0] = dat[1] + dat[9];
	buf[1] = dat[13] + dat[19];
	buf[2] = dat[39] + dat[47];
	buf[3] = dat[29] + dat[35];
	tmm0 = _mm_load_ps(buf);

	/* 5-8 */
	buf[0] = dat[7] + dat[15];
	buf[1] = dat[5] + dat[11];
	buf[2] = dat[33] + dat[41];
	buf[3] = dat[37] + dat[43];
	tmm1 = _mm_load_ps(buf);

	/* 9-12 */
	buf[0] = dat[0] + dat[8] + dat[16];
	buf[1] = dat[6] + dat[12] + dat[18];
	buf[2] = dat[32] + dat[40] + dat[48];
	buf[3] = dat[30] + dat[36] + dat[42];
	tmm2 = _mm_load_ps(buf);

	/* 13-16 */
	buf[0] = dat[2] + dat[3] + dat[4] + dat[10];
	buf[1] = dat[20] + dat[26] + dat[27] + dat[34];
	buf[2] = dat[38] + dat[44] + dat[45] + dat[46];
	buf[3] = dat[14] + dat[21] + dat[22] + dat[28];
	tmm3 = _mm_load_ps(buf);

	/* 16 13 14 15 */
	buf[4] = buf[3];
	buf[5] = buf[0];
	buf[6] = buf[1];
	buf[7] = buf[2];
	tmm6 = _mm_load_ps(&(buf[4]));

	/* 17-20 */
	buf[0] = dat[17];
	buf[1] = dat[25];
	buf[2] = dat[31];
	buf[3] = dat[23];
	tmm4 = _mm_load_ps(buf);

	/* 20 17 18 19 */
	buf[4] = buf[3];
	buf[5] = buf[0];
	buf[6] = buf[1];
	buf[7] = buf[2];
	tmm5 = _mm_load_ps(&(buf[4]));

	tmm4 = _mm_add_ps(tmm4, tmm0);
	tmm0 = _mm_add_ps(tmm0, tmm1);
	tmm0 = _mm_add_ps(tmm0, tmm2);
	tmm5 = _mm_add_ps(tmm5, tmm1);


	/*    E-H:tmm4
		  A-D:tmm0
		  I-L:tmm5
		 9-12:tmm2
		13-16:tmm3
		16 13 14 15:tmm6
		*/
	tmm0 = _mm_add_ps(tmm0, tmm3);
	tmm0 = _mm_add_ps(tmm0, tmm6);
	/* AA-DD:tmm0 */

	u.v = tmm0;
	buf[0] = u.a[3];
	buf[1] = u.a[0];
	buf[2] = u.a[1];
	buf[3] = u.a[2];
	tmm1 = _mm_load_ps(buf);
	/* DD AA BB CC:tmm1 */

	u.v = tmm4;
	buf[0] = u.a[2];
	buf[1] = u.a[3];
	buf[2] = u.a[0];
	buf[3] = u.a[1];
	tmm7 = _mm_load_ps(buf);
	/* G H E F: tmm7 */

	buf[0] = u.a[3];
	buf[1] = u.a[0];
	buf[2] = u.a[1];
	buf[3] = u.a[2];
	tmm4 = _mm_load_ps(buf);
	/* H E F G: tmm4 */

	u.v = tmm2;
	buf[0] = u.a[3];
	buf[1] = u.a[0];
	buf[2] = u.a[1];
	buf[3] = u.a[2];
	tmm3 = _mm_load_ps(buf);
	/* 12 9 10 11: tmm3 */
	mmm = buf[1] + buf[3] + dat[24];
	nnn = buf[2] + buf[0] + dat[24];

	/*		AA-DD: tmm0 
	  DD AA BB CC: tmm1
			 9-12:tmm2
	   12 9 10 11: tmm3
		  H E F G: tmm4
			  I-L:tmm5
		  G H E F: tmm7 */

	printf("AA-DD\n");
	printXmm(tmm0);
	printf("9-12\n");
	printXmm(tmm2);
	printf("I-L\n");
	printXmm(tmm5);
	printf("G H E F\n");
	printXmm(tmm7);

	tmm2 = _mm_add_ps(tmm2, tmm1);
	tmm2 = _mm_add_ps(tmm2, tmm5);
	/* DDD AAA BBB CCC: tmm2 */
	u.v = tmm2;
	aaa = u.a[1];
	bbb = u.a[2];
	ccc = u.a[3];
	ddd = u.a[0];

	tmm7 = _mm_add_ps(tmm7, tmm1);
	tmm7 = _mm_add_ps(tmm7, tmm5);
	/* HHH EEE FFF GGG: tmm7 */
	u.v = tmm7;
	eee = u.a[1];
	fff = u.a[2];
	ggg = u.a[3];
	hhh = u.a[0];

	tmm0 = _mm_add_ps(tmm0, tmm4);
	tmm0 = _mm_add_ps(tmm0, tmm3);
	/* III JJJ KKK LLL: tmm0 */
	u.v = tmm0;
	iii = u.a[0];
	jjj = u.a[1];
	kkk = u.a[2];
	lll = u.a[3];


	/*  DDD AAA BBB CCC: tmm3
		HHH EEE FFF GGG: tmm7
		III JJJ KKK LLL: tmm0
	*/

	u.v = tmm4;
	ooo = u.a[1] + u.a[3] + dat[24];
	ppp = u.a[0] + u.a[2] + dat[24];

	u.v = tmm5;
	qqq = u.a[0] + u.a[2] + dat[24];
	rrr = u.a[1] + u.a[3] + dat[24];

	/* ----------------stage2--------------- */
	float a5, b5, c5, d5, e5, f5, g5, h5, i5, j5, k5, l5,
		m5, n5, o5, p5, q5, r5;
	float c_21, c_22;
	c_21 = dat[23] + dat[24] + dat[25];
	c_22 = dat[17] + dat[24] + dat[31];

	/* 1-4 */
	buf[0] = dat[0] + dat[1] + dat[7] + dat[8];
	buf[1] = dat[5] + dat[6] + dat[12] + dat[13];
	buf[2] = dat[40] + dat[41] + dat[47] + dat[48];
	buf[3] = dat[35] + dat[36] + dat[42] + dat[43];
	tmm0 = _mm_load_ps(buf);

	buf[4] = buf[1];
	buf[5] = buf[2];
	buf[6] = buf[3];
	buf[7] = buf[0];
	tmm5 = _mm_load_ps(&(buf[4]));

	tmm0 = _mm_add_ps(tmm0, tmm5);

	/*  AA-DD: tmm0 */

	/* 5-8 */
	buf[0] = dat[3] + dat[10];
	buf[1] = dat[26] + dat[27];
	buf[2] = dat[38] + dat[45];
	buf[3] = dat[21] + dat[22];
	tmm1 = _mm_load_ps(buf);

	tmm0 = _mm_add_ps(tmm0, tmm1);

	buf[4] = buf[1];
	buf[5] = buf[2];
	buf[6] = buf[3];
	buf[7] = buf[0];
	tmm5 = _mm_load_ps(&(buf[4]));

	buf[4] = buf[3];
	buf[5] = buf[0];
	buf[6] = buf[1];
	buf[7] = buf[2];
	tmm6 = _mm_load_ps(&(buf[4]));

	q5 = buf[0] + buf[2] + c_22;
	r5 = buf[1] + buf[3] + c_21;


	/*  5-8: tmm1
	6 7 8 5: tmm5
	8 5 6 7: tmm6
		*/

	/* 9-12 */
	buf[0] = dat[2] + dat[9];
	buf[1] = dat[19] + dat[20];
	buf[2] = dat[39] + dat[46];
	buf[3] = dat[28] + dat[29];
	tmm2 = _mm_load_ps(buf);

	buf[4] = buf[1];
	buf[5] = buf[2];
	buf[6] = buf[3];
	buf[7] = buf[0];
	tmm7 = _mm_load_ps(&(buf[4]));

	m5 = buf[0] + buf[2] + c_22;
	n5 = buf[1] + buf[3] + c_22;
	/*     9-12: tmm2
	 10 11 12 9: tmm7
		*/

	/* 13-16 */
	buf[0] = dat[4] + dat[11];
	buf[1] = dat[33] + dat[34];
	buf[2] = dat[37] + dat[44];
	buf[3] = dat[14] + dat[15];
	tmm3 = _mm_load_ps(buf);

	buf[4] = buf[3];
	buf[5] = buf[0];
	buf[6] = buf[1];
	buf[7] = buf[2];
	tmm1 = _mm_load_ps(&(buf[4]));

	o5 = buf[0] + buf[2] + c_21;
	p5 = buf[1] + buf[3] + c_21;
	/*     13-16: tmm3
	 16 13 14 15: tmm1
	 */

	/* 17-20 */
	buf[0] = dat[16] + dat[17] + dat[18];
	buf[1] = dat[18] + dat[25] + dat[32];
	buf[2] = dat[30] + dat[31] + dat[32];
	buf[3] = dat[16] + dat[23] + dat[30];
	tmm4 = _mm_load_ps(buf);

	tmm4 = _mm_add_ps(tmm4, tmm0);
	tmm4 = _mm_add_ps(tmm4, tmm2);
	tmm4 = _mm_add_ps(tmm4, tmm3);

	/*    
	      A4-D4: tmm4
		6 7 8 5: tmm5
		8 5 6 7: tmm6
		   9-12: tmm2
	 10 11 12 9: tmm7
		  13-16: tmm3
	16 13 14 15: tmm1
		 */
	tmm6 = _mm_add_ps(tmm6, tmm4);
	tmm6 = _mm_add_ps(tmm6, tmm1);

	tmm5 = _mm_add_ps(tmm5, tmm4);
	tmm5 = _mm_add_ps(tmm5, tmm7);

	tmm1 = _mm_add_ps(tmm1, tmm4);
	tmm1 = _mm_add_ps(tmm1, tmm7);
	/*
		A5-D5: tmm6 
		E5-H5: tmm1
		I5-L5: tmm5
		*/
	u.v = tmm6;
	a5 = u.a[0];
	b5 = u.a[1];
	c5 = u.a[2];
	d5 = u.a[3];

	u.v = tmm1;
	e5 = u.a[0];
	f5 = u.a[1];
	g5 = u.a[2];
	h5 = u.a[3];

	u.v = tmm5;
	i5 = u.a[0];
	j5 = u.a[1];
	k5 = u.a[2];
	l5 = u.a[3];

	printf("group2\n");
	printf("a:%d c:%d n:%d\n", (int)a5, (int)c5, (int)n5);
	printf("b:%d d:%d m:%d\n", (int)b5, (int)d5, (int)m5);

	printf("e:%d g:%d r:%d\n", (int)e5, (int)g5, (int)r5);
	printf("f:%d h:%d q:%d\n", (int)f5, (int)h5, (int)q5);

	printf("i:%d k:%d p:%d\n", (int)i5, (int)k5, (int)p5);
	printf("j:%d l:%d o:%d\n", (int)j5, (int)l5, (int)o5);

	printf("group1\n");
	printf("a:%d c:%d p:%d\n", (int)aaa, (int)ccc, (int)ppp);
	printf("b:%d d:%d o:%d\n", (int)bbb, (int)ddd, (int)ooo);

	printf("e:%d g:%d n:%d\n", (int)eee, (int)ggg, (int)nnn);
	printf("f:%d h:%d m:%d\n", (int)fff, (int)hhh, (int)mmm);

	printf("i:%d k:%d r:%d\n", (int)iii, (int)kkk, (int)rrr);
	printf("j:%d l:%d q:%d\n", (int)jjj, (int)lll, (int)qqq);
		
	printDat(dat, width, height);
#if 0
	printRefLineVertical(dat, width, height);
#endif
	
}

void printDat(float* dat, int width, int height) {
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			printf("%02d ", (int)dat[i * width + j]);
		}
		printf("\n");
	}
	return;
}

void printRefLineVertical(float* dat, int width, int height) {
	float r, g, b;
	r = g = b = 0;

	for (int i = 0; i < 21; i++) g += dat[i];
	for (int i = 21; i < 28; i++) b += dat[i];
	for (int i = 28; i < 49; i++) r += dat[i];

	printf("g:%d\n b:%d\n r:%d\n", (int)g, (int)b, (int)r);
}