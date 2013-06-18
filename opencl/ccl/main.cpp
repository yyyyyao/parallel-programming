#include <math.h>
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <list>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>
#ifdef __APPLE__
#include <OpenCl/opencl.h>
#else
#include <CL/cl.h>
#endif

using namespace std;

list<unsigned int> getConnectPix8Dir(unsigned int cur, unsigned int rows, unsigned int cols, unsigned int step)
{
  list<unsigned int> l;
  if(cur % step != 0) {
    l.push_back(cur - 1);
    if(cur / step != 0)
      l.push_back(cur - 1 - step);
    if(cur / step != rows - 1)
      l.push_back(cur - 1 + step);
  }
  if(cur % step != cols - 1) {
    l.push_back(cur + 1);
    if(cur / step != 0)
      l.push_back(cur + 1 - step);
    if(cur / step != rows - 1)
      l.push_back(cur + 1 + step);
  }
  if(cur / step != 0)
    l.push_back(cur - step);
  if(cur / step != rows - 1)
    l.push_back(cur + step);

#if 0
  list<unsigned int>::iterator it = l.begin();
  while (it != l.end()) {
    ++it;
  }
#endif
  return l;
}

list<unsigned int> getConnectPix4Dir(unsigned int cur, unsigned int rows, unsigned int cols, unsigned int step)
{
  list<unsigned int> l;
  if(cur % step != 0)
    l.push_back(cur - 1);
  if(cur % step != cols - 1)
    l.push_back(cur + 1);
  if(cur / step != 0)
    l.push_back(cur - step);
  if(cur / step != rows - 1)
    l.push_back(cur + step);
  return l;
}

cv::Mat Ccl(const cv::Mat input, unsigned char threshold)
{
  int div  = 31;
  cv::Mat dir = input.clone();
  unsigned int cur, next;
  unsigned int cNum = 0;
  unsigned long long* index_mat;
  unsigned int index;
  unsigned int i, j;
  list<unsigned int> conList, tempList;
  index_mat = (unsigned long long*)malloc(dir.rows * dir.step * sizeof(unsigned long long));
  for(i = 0; i < dir.rows; i++) {
    for(j = 0; j < dir.cols; j++) {
      index = j + i * dir.step;
      index_mat[index] = index;
    }
  }

  for(i = 0; i < dir.rows; i++) {
    for(j = 0; j < dir.cols; j++) {
      index = j + i * dir.step;
      if (index_mat[index] < index)
        continue;
      conList.clear();
      conList.push_front(index);
      dir.data[index] = cNum * div % 255;
      while (!conList.empty()) {
        cur = conList.back();
        conList.pop_back();

        //tempList = getConnectPix8Dir(cur, input.rows, input.cols, input.step);
        tempList = getConnectPix4Dir(cur, input.rows, input.cols, input.step);
        while(!tempList.empty()) {
          next = tempList.back();
          tempList.pop_back();
          if (index_mat[next] <= index)
            continue;
          list<unsigned int>::iterator pos;
          pos = find(conList.begin(), conList.end(), next);
          if(pos != conList.end())
            continue;
          if(max(input.data[cur], input.data[next]) -
              min(input.data[cur], input.data[next]) > threshold)
            continue;

          dir.data[next] = cNum * div % 255;
          index_mat[next] = index;
          conList.push_front(next);
        }
      }
      cNum++;
    }
  }
  
  return dir;
}

bool compareMat(const cv::Mat mat1, const cv::Mat mat2)
{
  int err = 0;
  if(mat1.cols!= mat2.cols || mat1.rows != mat2.rows) return false;

  for(int i = 0; i < mat1.rows; i++) {
    for(int j = 0; j < mat1.cols; j++) {
      if(abs(mat1.data[j + i * mat1.step] - mat2.data[j + i * mat1.step]) > 1) {
        /*
        cout << "diffrent image!" << "cols:" << j << " rows:" << i << endl;
        cout << "value mat1:" << (int)mat1.data[j + i * mat1.step] << 
            " value mat2:" << (int)mat2.data[j + i * mat1.step] << endl;
        cout << "it's uncorrect!" << endl;
        */
        //return false;
        if(err++ > 10) return false;
      }
    }
  }
  cout << "it's correct!" << endl;
  return true;
}

int main(int argc, char* argv[])
{
  cv::Mat input;
  cv::Mat output;
  int threshold = 128;
 
  if(argc > 1) {
    input = cv::imread(argv[1], 0);
  }
  IplImage ipl_input = input;

  IplImage ipl_out = output;

  cv::Mat ans;
  ans = Ccl(input, threshold);

  //cvShowImage("", &ipl_out);
  cv::imwrite("./ans.bmp", ans);
  cv::imshow("", ans);
  cv::waitKey(0);
  return 0;
}
