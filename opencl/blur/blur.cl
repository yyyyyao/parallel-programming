__kernel void blur_cl(
  __global const uchar* src,
  __global uchar* dst,
  const int srcHeight,
  const int srcWidth,
  const int srcWidthStep,
  const int radius,
  const float fweight,
  __global float* f_dst)
{
  int i, j, k, l, _k, _l;
  float sum;

  i = get_global_id(0) / srcWidthStep;
  j = get_global_id(0) % srcWidthStep;

  sum = (j + i * srcWidthStep) % 255 / 3;
  dst[j + i * srcWidthStep] = convert_uchar_sat_rte(sum);
  f_dst[j + i * srcWidthStep] = fweight;
  return;
  //dst[j + i * srcWidthStep] = convert_uchar_sat_rte(sum);
  //dst[j + i * srcWidthStep] = src[j + i * srcWidthStep];

  //if(i >= srcHeight || j >= srcWidth) return;

  //i = get_global_id(0); /* rows */
  //j = get_global_id(1); /* cols */

  sum = 0;
  for(k = i - radius; k <= i + radius; k++) {
    if(k < 0) _k = k * -1;
    else if(k >= srcHeight) _k = (srcHeight- 1) * 2 - k;
    else _k = k;
    for(l = j - radius; l <= j + radius; l++) {
      if(l < 0) _l = l * -1;
      else if(l >= srcWidth) _l = (srcWidth - 1) * 2 - l;
      else _l = l;

      sum += src[_l + _k * srcWidthStep] * fweight;
    }
  }
  dst[j + i * srcWidthStep] = convert_uchar_sat_rte(sum);
  f_dst[j + i * srcWidthStep] = sum;
  return;
}

__kernel void __blur(
  __global uchar* src,
  __global uchar* dst,
  const int srcHeight,
  const int srcWidth,
  const int srcWidthStep,
  const int radius,
  const float fweight)
{
  int i;
  i = get_global_id(0);
  dst[i] = 111;
}
