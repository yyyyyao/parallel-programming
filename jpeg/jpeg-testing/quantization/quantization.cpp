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

void calcQuantTableByQuality(int* pSrc, int* pDst, int quality);
void applyQuantTable(const int* pSrc, int* pDst, int* quantTable);
bool compareMat(const cv::Mat mat1, const cv::Mat mat2);

static int s_std_lum_quant[64] = { 16,11,12,14,12,10,16,14,13,14,18,17,16,19,24,40,26,24,22,22,24,49,35,37,29,40,58,51,61,60,57,51,56,55,64,72,92,78,64,68,87,69,55,56,80,109,81,87,95,98,103,104,103,62,77,113,121,112,100,120,92,101,103,99 };
static int s_std_croma_quant[64] = { 17,18,18,24,21,24,47,26,26,47,99,66,56,66,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99 };

static const int dctMcu[64] = {
223, 130, 40, 16, 11, 8, 2, -1, 
41, -34, -14, -10, -4, 0, -1, 3, 
-7, 10, -12, 2, 2, -5, 1, -1, 
22, -7, 9, 2, 0, 1, -3, 2, 
-7, 4, -6, 3, -1, -2, 4, -1, 
5, 2, -1, -4, 0, 1, -1, -1, 
4, -5, 3, -1, 0, 2, 0, -1, 
-5, 5, -2, 3, 0, -2, 1, -1
};

int main(int argc, char* argv[])
{
  int quantTable[64];
  int result[64];
  memset(quantTable, 0, sizeof(int) * 64);
  memset(result, 0, sizeof(int) * 64);

  printf("-----dct value-----\n");
  for(int i = 0; i < 8; i++) {
    for (int j = 0; j < 8; j++) {
      printf("%4d ", dctMcu[i * 8 + j]);
    }
    printf("\n");
  }
  calcQuantTableByQuality(s_std_lum_quant, quantTable, 80);
  applyQuantTable(dctMcu, result, quantTable);

  printf("------quantization table------\n");
  for(int i = 0; i < 8; i++) {
    for (int j = 0; j < 8; j++) {
      printf("%2d ", quantTable[i * 8 + j]);
    }
    printf("\n");
  }

  printf("------quqntized data-----\n");
  for(int i = 0; i < 8; i++) {
    for (int j = 0; j < 8; j++) {
      printf("%2d ", result[i * 8 + j]);
    }
    printf("\n");
  }
  return 0;
}

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
