/*
 * =====================================================================================
 *
 *       Filename:  blur.cpp
 *
 *     
 *
 *        Version:  1.0
 *        Created:  2013/07/28 14時20分48秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tetsuhiko Yao (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <iostream>
#include <string>
#include <stdio.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#define uint8 unsigned char
#define uint16 unsigned short
#define uint unsigned int
using namespace std;

static void RGB_to_YCC(uint8* , const uint8 *, int );

int main(int argc, char* argv[])
{
  const string windowName = "result window";
  cv::Mat input;
  cv::Mat refInput;
  cv::Size ksize;
 
  if(argc > 1) {
    input = cv::imread(argv[1], 1);
    refInput = cv::imread(argv[1], 0);
  }
  int height = input.rows;
  int width = input.cols;
  int comp = input.elemSize();

  printf("height:%d width:%d comp:%d\n", height, width, comp);

  cv::Mat output(input.rows, input.cols, CV_8UC3);
  cv::Mat grayOutput(input.rows, input.cols, CV_8UC1);

  ksize.width = 3;
  ksize.height = 3;
  if(argc > 2) {
    ksize.width = atoi(argv[2]);
    ksize.height = atoi(argv[2]);
  }

  unsigned char* srcPtr = (unsigned char*)input.data;
  unsigned char* refPtr = (unsigned char*)refInput.data;
  unsigned char* dstPtr = (unsigned char*)output.data;
  unsigned char* gDstPtr = (unsigned char*)grayOutput.data;
  unsigned char* srcPtrBuf;
  unsigned char* dstPtrBuf;

  /* execute refrence codes */

  for(int i = 0; i < height; i++) {
    srcPtrBuf = srcPtr + i * width * comp;
    dstPtrBuf = dstPtr + i * width * comp;
    RGB_to_YCC(dstPtrBuf, srcPtrBuf, width);
  }

  /* check Y value. */
  for(int i = 0; i < height; i++) {
    for(int j = 0; j < width; j++) {
      gDstPtr[i * width + j] = dstPtr[i * width * comp + j * comp];
    }
  }

  cv::namedWindow(windowName, CV_WINDOW_AUTOSIZE);
  cv::imshow(windowName, grayOutput);
  cv::waitKey(0);

  return 0;
}

const int YR = 19595, YG = 38470, YB = 7471, CB_R = -11059, CB_G = -21709, CB_B = 32768, CR_R = 32768, CR_G = -27439, CR_B = -5329;
static inline uint8 clamp(int i) { if (static_cast<uint>(i) > 255U) { if (i < 0) i = 0; else if (i > 255) i = 255; } return static_cast<uint8>(i); }

static void RGB_to_YCC(uint8* pDst, const uint8 *pSrc, int num_pixels)
{
  for ( ; num_pixels; pDst += 3, pSrc += 3, num_pixels--)
  {
    const int r = pSrc[0], g = pSrc[1], b = pSrc[2];
    pDst[0] = static_cast<uint8>((r * YR + g * YG + b * YB + 32768) >> 16);
    pDst[1] = clamp(128 + ((r * CB_R + g * CB_G + b * CB_B + 32768) >> 16));
    pDst[2] = clamp(128 + ((r * CR_R + g * CR_G + b * CR_B + 32768) >> 16));
  }
}
