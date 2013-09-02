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
#define int8 char
#define uint16 unsigned short
#define uint unsigned int
using namespace std;

  const int8 inputMcu[64] ={
  57,49,44,39,33,28,23,22,
  55,45,39,33,28,21,21,19,
  55,44,37,31,20,17,17,16,
  58,44,37,29,17,15,16,14,
  66,46,31,22,16,14,12,10,
  71,49,32,24,14,12,21,19,
  69,50,30,22,12,4,-1,-2,
  62,42,25,17,7,1,-16,-16
  };
  int output[64]; 

static void twoDimIdct(const int* pSrc, int* pDst);
static void twoDimDct(const int8* pSrc, int* pDst);
static void refTwoDimDct(const int8* pSrc, int* pDst);
bool compareMat(const cv::Mat mat1, const cv::Mat mat2);

int main(int argc, char* argv[])
{
  refTwoDimDct(inputMcu, output);
  for(int i = 0; i < 8; i++) {
    for(int j = 0; j < 8; j++) {
      printf("%d, ", output[i * 8 + j]);
    }
    printf("\n");
  }
  memset(output, 0, sizeof(int) * 64);
  printf("--------------------------\n");
  twoDimDct(inputMcu, output);

  for(int i = 0; i < 8; i++) {
    for(int j = 0; j < 8; j++) {
      printf("%5d ", output[i * 8 + j]);
    }
    printf("\n");
  }

  int result[64];
  twoDimIdct(output, result);

  printf("--------------------------\n");
  for(int i = 0; i < 8; i++) {
    for(int j = 0; j < 8; j++) {
      printf("%5d ", result[i * 8 + j]);
    }
    printf("\n");
  }

  printf("--------------------------\n");
  for(int i = 0; i < 8; i++) {
    for(int j = 0; j < 8; j++) {
      printf("%5d ", inputMcu[i * 8 + j]);
    }
    printf("\n");
  }
  return 0;
}

void twoDimIdct(const int* pSrc, int* pDst)
{
  double cu, cv;
  double uCof, vCof;
  double cucv, val;

  for(int y = 0; y < 8; y++) {
    for(int x = 0; x < 8; x++) {
      val = 0;
      for(int v = 0; v < 8; v++) {
        if(v) cv = 1;
        else cv = 1 / sqrt(2);
        vCof = (2 * y + 1) * v * M_PI / (2 * 8);
        for(int u = 0; u < 8; u++) {
          if(u) cu = 1;
          else cu = 1 / sqrt(2);
          uCof = (2 * x + 1) * u * M_PI / (2 * 8);
          val += cu * cv * pSrc[v * 8 + u] * cos(uCof) * cos(vCof);
        }
      }
      pDst[y * 8 + x] = floor((2 * val / 8) + 0.5);
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

void refTwoDimDct(const int8* pSrc, int* pDst)
{
  double cv, cu;
  double sum;
  int after[8][8];
  for(int v =0; v < 8; v++)
  {
    if(v) cv =1.0;
    if(!v) cv = 1 / sqrt(2.0);
    for(int u = 0; u < 8; u++)
    {
      if(u) cu =1.0;
      if(!u) cu = 1 / sqrt(2.0);
      sum = 0.0;
      for(int y = 0; y < 8; y++)
      {
        for(int x = 0; x < 8; x++)
        {
          double a = cos((2 * x + 1)*u*M_PI/16) * cos(( 2 * y + 1)*v*M_PI/16);
          double d = pSrc[y * 8 + x];
          sum += d * a;
        }
      }
      pDst[v * 8 + u] = int(floor(cu*cv*sum/4 + 0.5));
    }
  }
  return;
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
