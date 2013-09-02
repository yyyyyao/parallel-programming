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

static void blockSplitMcuLine(const uint8* pSrc, uint8* pDst, int width);
static void blockSplitting(const uint8* pSrc, uint8* pDst, int width, int height);
static void reverseBlockSplitting(const uint8* pSrc, uint8* pDst, int width, int height);
bool compareMat(const cv::Mat mat1, const cv::Mat mat2);

int main(int argc, char* argv[])
{
  const string windowName = "result window";
  cv::Mat input;
  cv::Mat refInput;
  cv::Size ksize;
 
  if(argc > 1) {
    input = cv::imread(argv[1], 0);
  }
  int height = input.rows;
  int width = input.cols;
  int comp = input.elemSize();

  printf("height:%d width:%d comp:%d\n", height, width, comp);

  cv::Mat output(input.rows, input.cols, CV_8UC1);
  cv::Mat ref(input.rows, input.cols, CV_8UC1);
  cv::Mat result(input.rows, input.cols, CV_8UC1);

  ksize.width = 3;
  ksize.height = 3;
  if(argc > 2) {
    ksize.width = atoi(argv[2]);
    ksize.height = atoi(argv[2]);
  }

  unsigned char* srcPtr = (unsigned char*)input.data;
  unsigned char* dstPtr = (unsigned char*)output.data;
  unsigned char* refPtr = (unsigned char*)ref.data;
  unsigned char* resultPtr = (unsigned char*)result.data;
  unsigned char* srcPtrBuf;
  unsigned char* dstPtrBuf;

  /* execute refrence codes */

  for(int i = 0; i < height; i+=8) {
    srcPtrBuf = srcPtr + i * width * comp;
    dstPtrBuf = dstPtr + i * width * comp;
    blockSplitMcuLine(srcPtrBuf, dstPtrBuf, width);
  }

  blockSplitting(srcPtr, refPtr, width, height);
  reverseBlockSplitting(dstPtr, resultPtr, width, height);

  compareMat(ref, output);
  compareMat(result, input);
  cv::namedWindow(windowName, CV_WINDOW_AUTOSIZE);
  cv::imshow(windowName, result);
  cv::waitKey(0);

  return 0;
}

static void blockSplitting(const uint8* pSrc, uint8* pDst, int width, int height)
{
  const int blockWidth = width / 8;
  for(int y = 0; y < height; y++) {
    for(int x = 0; x < width; x++) {
      int dstIndex = ((x % width) / 8 + (y / 8) * blockWidth) * 64 +
        (y % 8) * 8 + x % 8;
      int srcIndex = x + y * width;
      pDst[dstIndex] = pSrc[srcIndex];
    }
  }
}

static void reverseBlockSplitting(const uint8* pSrc, uint8* pDst, int width, int height)
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

static void blockSplitMcuLine(const uint8* pSrc, uint8* pDst, int width)
{
  for(int i = 0; i < 8; i++) { /* because MUC height is 8. */
    for(int j = 0; j < width; j++) {
      int insideMcuIndex = i * 8 + j % 8;
      int mcuNum = j / 8;
      pDst[mcuNum * 64 + insideMcuIndex] = pSrc[i * width + j];
    }
  }
}

bool compareMat(const cv::Mat mat1, const cv::Mat mat2)
{
  int err = 0;
  if(mat1.cols!= mat2.cols || mat1.rows != mat2.rows) return false;

  for(int i = 0; i < mat1.rows; i++) {
    for(int j = 0; j < mat1.cols; j++) {
      if(abs(mat1.data[j + i * mat1.step] - mat2.data[j + i * mat1.step]) > 1) {
        cout << "diffrent image!" << "cols:" << j << " rows:" << i << endl;
        cout << "value mat1:" << (int)mat1.data[j + i * mat1.step] << 
            " value mat2:" << (int)mat2.data[j + i * mat1.step] << endl;
        cout << "it's uncorrect!" << endl;
        //return false;
        if(err++ > 10) return false;
      }
    }
  }
  cout << "it's correct!" << endl;
  return true;
}
