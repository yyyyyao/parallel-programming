/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2013年09月11日 18時46分41秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <stdio.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "jpgRef.hpp"

using namespace std;
int main(int argc, char* argv[])
{
  cv::Mat input;
  cv::Mat refInput;
  int quality = 90;
  char* outputFile = (char*)"result.jpg";
  if(argc > 1) {
    input = cv::imread(argv[1], 1);
    refInput = cv::imread(argv[1], 0);
  }

  if(argc > 2)
    outputFile = argv[2];
  if(argc > 3)
    quality = atoi(argv[3]);

  jpgRef j(input, outputFile, quality);
  j.execute();

  //cv::namedWindow("ref", CV_WINDOW_AUTOSIZE);
  //cv::imshow("ref", refInput);
  //cv::waitKey(0);
  return 0;
}
