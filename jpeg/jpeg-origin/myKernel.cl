//#define NO_SUBSAMPLING 1
__kernel void colConv_function(__global unsigned char* src,
    __global unsigned char* dst,
    const int imgHeight,
    const int imgWidth,
    __global unsigned char* dstY,
    __global unsigned char* dstCb,
    __global unsigned char* dstCr
    ) {
  const int id0 = get_global_id(0);
  const int id1 = get_global_id(1);
  /* one thread executes one pixel. */
  if(id0 > imgHeight) return;
  if(id1 > imgWidth) return;

  int i = id0 * imgWidth + id1;

  const int imgSize = imgWidth * imgHeight;

  float b = src[i];
  float g = src[imgSize + i];
  float r = src[2*imgSize + i];

  float y  =   0 + 0.299 * r + 0.587 * g + 0.114 * b;
  float cb = 128 - 0.169 * r - 0.331 * g + 0.500 * b;
  float cr = 128 + 0.500 * r - 0.419 * g - 0.081 * b;

  dstY[i] = (unsigned char)y;
#ifdef NO_SUBSAMPLING
  dstCb[i] = (unsigned char)cb;
  dstCr[i] = (unsigned char)cr;

#else
  const int halfImgWidth = imgWidth / 2;
  if(id0 % 2 == 0 && id1 % 2 == 0) {
    int subIndex = (id0 / 2) * halfImgWidth + (id1 / 2);
    dstCb[subIndex] = (unsigned char)cb;
    dstCr[subIndex] = (unsigned char)cr;
  }
#endif
}
