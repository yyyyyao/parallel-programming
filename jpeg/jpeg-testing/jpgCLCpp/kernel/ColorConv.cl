__kernel void ColorConv_SubSampling(
    global unsigned char* src,
    const int imgHeight,
    const int imgWidth,
    __global unsigned char* dstY,
    __global unsigned char* dstCb,
    __global unsigned char* dstCr) {
  int gId0 = get_global_id(0); /* width / 2 */
  int gId1 = get_global_id(1); /* height / 2 */
  int i, j;
  float r, g, b, y, cb, cr;

  int baseIndex = (gId1 * 2) * imgWidth * 3 + (gId0 * 2 * 3);

  for(i = 0; i < 2; i++) {
    for(j = 0; j < 2; j++) {
      r = src[baseIndex + i * imgWidth * 3 + j * 3 + 2];
      g = src[baseIndex + i * imgWidth * 3 + j * 3 + 1];
      b = src[baseIndex + i * imgWidth * 3 + j * 3];
      y = 0.299f * r + 0.587f * g + 0.114f * b;
      dstY[(gId1 * 2 + i) * imgWidth + gId0 * 2 + j] = (unsigned char)y;

      if(i == 0 && j == 0) { /* top-left pixel */
        cb = 128.0f - 0.169f * r - 0.331f * g + 0.500f * b;
        cr = 128.0f + 0.500f * r - 0.419f * g - 0.081f * b;
        dstCb[gId1 * (imgWidth >> 1) + gId0] = cb;
        dstCr[gId1 * (imgWidth >> 1) + gId0] = cr;
      }
    }
  }
}
