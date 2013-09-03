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

#define int8 char
#define uint8 unsigned char
#define uint16 unsigned short
#define uint unsigned int
using namespace std;

static int s_std_lum_quant[64] = { 16,11,12,14,12,10,16,14,13,14,18,17,16,19,24,40,26,24,22,22,24,49,35,37,29,40,58,51,61,60,57,51,56,55,64,72,92,78,64,68,87,69,55,56,80,109,81,87,95,98,103,104,103,62,77,113,121,112,100,120,92,101,103,99 };
static int s_std_croma_quant[64] = { 17,18,18,24,21,24,47,26,26,47,99,66,56,66,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99 };

const int YR = 19595, YG = 38470, YB = 7471, CB_R = -11059, CB_G = -21709, CB_B = 32768, CR_R = 32768, CR_G = -27439, CR_B = -5329;
static inline uint8 clamp(int i) { if (static_cast<uint>(i) > 255U) { if (i < 0) i = 0; else if (i > 255) i = 255; } return static_cast<uint8>(i); }

void applyQuantTable(const int* pSrc, int* pDst, int* quantTable)
{
  for(int i = 0; i < 8; i++) {
    for(int j = 0; j < 8; j++) {
      pDst[i * 8 + j] = pSrc[i * 8 + j] / quantTable[i * 8 + j];
    }
  }
}

void calcQuantTableByQuality(int* pSrc, int* pDst, int quality)
{
    int s = (quality < 50) ? (5000 / quality) : (200 - (2 * quality));
    for ( int i = 0; i < 64; i++ ) {
        int value = (s * pSrc[i] + 50) / 100;
        if ( value == 0 ) {
            value = 1;
        }
        if ( value > 255 ) {
            value = 255;
        }
        pDst[i] = value;
    }
}
void centerization(const uint8* pSrc, int8* pDst, int width, int height)
{
  for(int i = 0; i < height; i++) {
    for(int j = 0; j < width; j++) {
      pDst[i * width + j] = (int8)(pSrc[i * width + j] - 128);
    }
  }

}
void twoDimDct(const int8* pSrc, int* pDst)
{
  double cu, cv;
  double xCof, yCof;
  double cucv, val;
  for(int v = 0;v < 8; v++) {
    if(v) cv = 1;
    else cv = 1 / sqrt(2);
    for(int u = 0; u < 8; u++) {
      val = 0;
      if(u) cu = 1;
      else cu = 1 / sqrt(2);
      for(int y = 0; y < 8; y++) {
        yCof = (2 * y + 1) * v * M_PI / (2 * 8);
        for(int x = 0; x < 8; x++) {
          xCof = (2 * x + 1) * u * M_PI / (2 * 8);
          val += (int)(pSrc[y * 8 + x]) * cos(xCof) * cos(yCof);
        }
      }
      cucv = (2 * cu * cv) / 8;
      pDst[v * 8 + u] = floor(val * cucv + 0.5);
    }
  }
}
void reverseBlockSplitting(const uint8* pSrc, uint8* pDst, int width, int height)
{
  const int blockWidth = width / 8;
  for(int y = 0; y < height; y++) {
    for(int x = 0; x < width; x++) {
      int dstIndex = ((x % width) / 8 + (y / 8) * blockWidth) * 64 +
        (y % 8) * 8 + x % 8;
      int srcIndex = x + y * width;
      pDst[srcIndex] = pSrc[dstIndex];
    }
  }
}
void blockSplitMcuLine(const uint8* pSrc, uint8* pDst, int width)
{
  for(int i = 0; i < 8; i++) { /* because MUC height is 8. */
    for(int j = 0; j < width; j++) {
      int insideMcuIndex = i * 8 + j % 8;
      int mcuNum = j / 8;
      pDst[mcuNum * 64 + insideMcuIndex] = pSrc[i * width + j];
    }
  }
}
void SubSamplingH2V2(const uint8* pSrc, uint8* pDst, int width, int height)
{
  int halfHeight = ceil(height / 2);
  int halfWidth = ceil(width / 2);

  for(int i = 0; i < halfHeight; i++) {
    for(int j = 0; j < halfWidth; j++) {
      pDst[i * halfWidth + j] = pSrc[(i * 2) * width + (j * 2)];
    }
  }
}
void deSplitFromComp(uint8* pSrc, const uint8* pY, const uint8* pCb, const uint8* pCr,
    int height, int width)
{
  for(int i = 0; i < height; i++) {
    for(int j = 0; j < width; j++) {
      pSrc[i * width * 3 + j * 3 + 0] = pY[i * width + j] ;
      pSrc[i * width * 3 + j * 3 + 1] = pCb[i * width + j];
      pSrc[i * width * 3 + j * 3 + 2] = pCr[i * width + j];
    }
  }
}

void splitToComp(const uint8* pSrc, uint8* pY, uint8* pCb, uint8* pCr, int height, int width)
{
  for(int i = 0; i < height; i++) {
    for(int j = 0; j < width; j++) {
      pY[i * width + j] = pSrc[i * width * 3 + j * 3 + 0];
      pCb[i * width + j] = pSrc[i * width * 3 + j * 3 + 1];
      pCr[i * width + j] = pSrc[i * width * 3 + j * 3 + 2];
    }
  }
}

static void RGB_to_YCC(uint8* pDst, const uint8 *pSrc, int num_pixels)
{
  for ( ; num_pixels; pDst += 3, pSrc += 3, num_pixels--)
  {
    const int r = pSrc[2], g = pSrc[1], b = pSrc[0];
    pDst[0] = static_cast<uint8>((r * YR + g * YG + b * YB + 32768) >> 16);
    pDst[1] = clamp(128 + ((r * CB_R + g * CB_G + b * CB_B + 32768) >> 16));
    pDst[2] = clamp(128 + ((r * CR_R + g * CR_G + b * CR_B + 32768) >> 16));
  }
}

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
  int halfWidth = ceil(width / 2);
  int halfHeight = ceil(height / 2);
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

  /*  RGB to YCbCr. */
  for(int i = 0; i < height; i++) {
    srcPtrBuf = srcPtr + i * width * comp;
    dstPtrBuf = dstPtr + i * width * comp;
    RGB_to_YCC(dstPtrBuf, srcPtrBuf, width);
  }

  /* split data for each comp. */
  cv::Mat Y(input.rows, input.cols, CV_8UC1);
  unsigned char* YPtr = (unsigned char*)Y.data;
  cv::Mat Cb(input.rows, input.cols, CV_8UC1);
  unsigned char* CbPtr = (unsigned char*)Cb.data;
  cv::Mat Cr(input.rows, input.cols, CV_8UC1);
  unsigned char* CrPtr = (unsigned char*)Cr.data;

  splitToComp(dstPtr, YPtr, CbPtr, CrPtr, width, height);

  /* SubSampling H2V2 */
  cv::Mat subCb(halfHeight, halfWidth, CV_8UC1);
  unsigned char* subCbPtr = (unsigned char*)subCb.data;
  cv::Mat subCr(halfHeight, halfWidth, CV_8UC1);
  unsigned char* subCrPtr = (unsigned char*)subCr.data;
  SubSamplingH2V2(CbPtr, subCbPtr, width, height);
  SubSamplingH2V2(CrPtr, subCrPtr, width, height);

  cv::Mat blockY(input.rows, input.cols, CV_8UC1);
  unsigned char* blockYPtr = (unsigned char*)blockY.data;
  unsigned char* YPtrBuf;
  unsigned char* blockYPtrBuf;
  /* block Splitting */
  for(int i = 0; i < height; i+=8) {
    YPtrBuf = YPtr + i * width ;
    blockYPtrBuf = blockYPtr + i * width;
    blockSplitMcuLine(YPtrBuf, blockYPtrBuf, width);
  }

  cv::Mat blockCb(halfHeight, halfWidth, CV_8UC1);
  unsigned char* blockCbPtr = (unsigned char*)blockCb.data;
  unsigned char* subCrPtrBuf;
  unsigned char* blockCrPtrBuf;
  cv::Mat blockCr(halfHeight, halfWidth, CV_8UC1);
  unsigned char* blockCrPtr = (unsigned char*)blockCr.data;
  unsigned char* subCbPtrBuf;
  unsigned char* blockCbPtrBuf;

  /* block Splitting */
  for(int i = 0; i < halfHeight; i+=8) {
    subCbPtrBuf = subCbPtr + i * halfWidth ;
    blockCbPtrBuf = blockCbPtr + i * halfWidth;
    blockSplitMcuLine(subCbPtrBuf, blockCbPtrBuf, halfWidth);

    subCrPtrBuf = subCrPtr + i * halfWidth ;
    blockCrPtrBuf = blockCrPtr + i * halfWidth;
    blockSplitMcuLine(subCrPtrBuf, blockCrPtrBuf, halfWidth);
  }

  /* apply reverse block splitting */
  cv::Mat reverseBlockY(input.rows, input.cols, CV_8UC1);
  unsigned char* reverseBlockYPtr = (unsigned char*)reverseBlockY.data;
  reverseBlockSplitting(blockYPtr, reverseBlockYPtr, width, height);

  cv::Mat reverseBlockCb(halfHeight, halfWidth, CV_8UC1);
  unsigned char* reverseBlockCbPtr = (unsigned char*)reverseBlockCb.data;
  reverseBlockSplitting(blockCbPtr, reverseBlockCbPtr, halfWidth, halfHeight);

  /* Centerization for DCT */
  cv::Mat centerY(input.rows, input.cols, CV_8UC1);
  char* centerYPtr= (char*)centerY.data;
  centerization(blockYPtr, centerYPtr, width, height);

  cv::Mat centerCb(halfHeight, halfWidth, CV_8UC1);
  char* centerCbPtr= (char*)centerCb.data;
  centerization(blockCbPtr, centerCbPtr, halfWidth, halfHeight);

  cv::Mat centerCr(halfHeight, halfWidth, CV_8UC1);
  char* centerCrPtr= (char*)centerCr.data;
  centerization(blockCrPtr, centerCrPtr, halfWidth, halfHeight);

  /* 2D DCT */
  cv::Mat dctY(input.rows, input.cols, CV_32SC1);
  int* dctYPtr= (int*)dctY.data;
  cv::Mat dctCr(halfHeight, halfWidth, CV_32SC1);
  int* dctCrPtr= (int*)dctCr.data;
  cv::Mat dctCb(halfHeight, halfWidth, CV_32SC1);
  int* dctCbPtr= (int*)dctCb.data;

  char* centerYPtrBuf;
  char* centerCbPtrBuf;
  char* centerCrPtrBuf;
  int* dctYPtrBuf;
  int* dctCbPtrBuf;
  int* dctCrPtrBuf;

  for(int i = 0; i < height; i+=8) {
    for(int j = 0; j < width; j+=8) {
      centerYPtrBuf = centerYPtr + i * width + j;
      dctYPtrBuf = dctYPtr + i * width + j;
      twoDimDct(centerYPtrBuf, dctYPtrBuf);
    }
  }

  for(int i = 0; i < halfHeight; i+=8) {
    for(int j = 0; j < halfWidth; j+=8) {
      centerCbPtrBuf = centerCbPtr + i * halfWidth + j;
      dctCbPtrBuf = dctCbPtr + i * halfWidth + j;
      twoDimDct(centerCbPtrBuf, dctCbPtrBuf);

      centerCrPtrBuf = centerCrPtr + i * halfWidth + j;
      dctCrPtrBuf = dctCrPtr + i * halfWidth + j;
      twoDimDct(centerCrPtrBuf, dctCrPtrBuf);
    }
  }

  /* Quantization */
  int quantTableY[64];
  int quantTableCbCr[64];
  int quality = 80;
  calcQuantTableByQuality(s_std_lum_quant, quantTableY, quality);
  calcQuantTableByQuality(s_std_croma_quant, quantTableCbCr, quality);

  cv::Mat quantY(input.rows, input.cols, CV_32SC1);
  int* quantYPtr= (int*)quantY.data;
  cv::Mat quantCr(halfHeight, halfWidth, CV_32SC1);
  int* quantCrPtr= (int*)quantCr.data;
  cv::Mat quantCb(halfHeight, halfWidth, CV_32SC1);
  int* quantCbPtr= (int*)quantCb.data;
  int* quantYPtrBuf;
  int* quantCrPtrBuf;
  int* quantCbPtrBuf;


  for(int i = 0; i < height; i+=8) {
    for(int j = 0; j < width; j+=8) {
      dctYPtrBuf = dctYPtr + i * width + j;
      quantYPtrBuf = quantYPtr + i * width + j;
      applyQuantTable(dctYPtrBuf, quantYPtrBuf, quantTableY);
    }
  }
  for(int i = 0; i < halfHeight; i+=8) {
    for(int j = 0; j < halfWidth; j+=8) {
      dctCbPtrBuf = dctCbPtr + i * halfWidth + j;
      quantCbPtrBuf = quantCbPtr + i * halfWidth + j;
      applyQuantTable(dctCbPtrBuf, quantCbPtrBuf, quantTableCbCr);

      dctCrPtrBuf = dctCrPtr + i * halfWidth + j;
      quantCrPtrBuf = quantCrPtr + i * halfWidth + j;
      applyQuantTable(dctCrPtrBuf, quantCrPtrBuf, quantTableCbCr);
    }
  }

  /* Entropy Encoding */

  /* Add Header */

  cv::namedWindow(windowName, CV_WINDOW_AUTOSIZE);
  cv::imshow(windowName, centerY);
  cv::waitKey(0);

  return 0;
}
