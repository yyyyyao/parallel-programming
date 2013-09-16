__kernel void transpose(
  __global const unsigned char* src,
  __global unsigned char* dst,
  const int srcHeight,
  const int srcWidth,
  const int srcWidthStep)
{
  int i = get_global_id(0);
  int j = get_global_id(1);

  if( i < srcHeight && j < srcWidth) {
    dst[get_global_id(0) * srcWidthStep + get_global_id(1)] = 
        src[get_global_id(1) * srcWidthStep + get_global_id(0)];
  }

}

__kernel void blur_cl_local(
  __global const uchar* src,
  __global uchar* dst,
  const int srcHeight,
  const int srcWidth,
  const int srcWidthStep,
  const int radius,
  const float fweight,
  __local uchar* localImg)

{
  int i, j, k, l, _k, _l;
  float sum = 0;

  i = get_global_id(0);
  j = get_global_id(1);

  int workItemStartRow = get_group_id(0) * get_local_size(0);
  int workItemStartCol = get_group_id(1) * get_local_size(1);

  int localRowIndex = get_local_id(0);
  int localColIndex = get_local_id(1);

  int localMemHeight = get_local_size(0) + radius * 2;
  int localMemWidth = get_local_size(1) + radius * 2;

  for(k = get_global_id(0) - radius;
      k < (signed int)(workItemStartRow + get_local_size(0) + radius);
      k += get_local_size(0)) {
    if(k < 0) _k = k * -1;
    else if(k >= srcHeight) _k = (srcHeight- 1) * 2 - k;
    else _k = k;

    for(l = get_global_id(1) - radius;
        l < (signed int)(workItemStartCol + get_local_size(1) + radius);
        l += get_local_size(1)) {
      if(l < 0) _l = l * -1;
      else if(l >= srcWidth) _l = (srcWidth - 1) * 2 - l;
      else _l = l;

      int localImgIndex = (k - workItemStartRow + radius) * localMemWidth +
          (l - workItemStartCol + radius);
      int globalImgIndex = _k * srcWidthStep + _l;
      localImg[localImgIndex] = src[globalImgIndex];
    }
  }

  barrier(CLK_LOCAL_MEM_FENCE);

  //dst[i * srcWidthStep + j] = 
  //    localImg[(localRowIndex + radius) * localMemWidth + localColIndex + radius];
  //return;

  for(k = localRowIndex; k <= localRowIndex + radius * 2; k++) {
    for(l = localColIndex; l <= localColIndex + radius * 2; l++) {
      sum += localImg[l + k * localMemWidth] * fweight;
    }
  }

  dst[j + i * srcWidthStep] = convert_uchar_sat_rte(sum);
  return;
}

__kernel void blur_cl_naive(
  __global const uchar* src,
  __global uchar* dst,
  const int srcHeight,
  const int srcWidth,
  const int srcWidthStep,
  const int radius,
  const float fweight)
{
  int i, j, k, l, _k, _l;
  float sum;

  i = get_global_id(0);
  j = get_global_id(1);

  //dst[get_global_id(1) * srcWidthStep + get_global_id(0)] = 
  //    src[get_global_id(1) * srcWidthStep + get_global_id(0)];
  //return;

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
  return;
}
