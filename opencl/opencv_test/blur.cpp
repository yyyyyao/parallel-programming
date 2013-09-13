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
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace std;

int main(int argc, char* argv[])
{
  const string windowName = "result window";
  cv::Mat input;
  cv::Mat output;
  cv::Size ksize;
 
  if(argc > 1) {
    input = cv::imread(argv[1], 0);
  }
  IplImage ipl = input;

  ksize.width = 3;
  ksize.height = 3;
  if(argc > 2) {
    ksize.width = atoi(argv[2]);
    ksize.height = atoi(argv[2]);
  }

  /* execute refrence codes */
  blur(input, output, ksize, cv::Point(-1, -1), cv::BORDER_REFLECT_101);

  cv::namedWindow(windowName, CV_WINDOW_AUTOSIZE);
  cv::imshow(windowName, output);
  cv::waitKey(0);

  return 0;
}
