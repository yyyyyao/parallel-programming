#include <math.h>
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <CL/cl.h>

#define MAX_SOURCE_SIZE (0x100000)

using namespace std;

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

void myblur(const cv::Mat& input, cv::Mat& dst, cv::Size size)
{
  int i, j, k, l, _k, _l;
  int radius;
  float sum, fweight;

  if(size.width % 2 == 0 ||
      size.width != size.height) {
    cout << "please set kernel.width/height to odd!" << endl;
    return;
  }

  radius = (size.width -1) / 2;
  fweight = (float)1 / (size.width * size.height);

  for(i = 0; i < input.rows; i++) {
    for(j = 0; j < input.cols; j++) {
      sum = 0;
      for(k = i - radius; k <= i + radius; k++) {
        if(k < 0) _k = k * -1;
        else if(k >= input.rows) _k = (input.rows - 1) * 2 - k;
        else _k = k;
        for(l = j - radius; l <= j + radius; l++) {
          if(l < 0) _l = l * -1;
          else if(l >= input.cols) _l = (input.cols - 1) * 2 - l;
          else _l = l;

          sum += input.data[_l + _k * input.step] * fweight;
#if 0
          if(i == 0 && j == 0) {
            cout << "sum:" << sum << endl;
            cout << "input.data:" << (int)src.data[_l + _k * src.step] << endl;
            cout << "l:" << l << " k:" << k << "_l:" << _l << "_k" << _k <<endl;
          }
#endif
        }
      }
      dst.data[j + i * input.step] = sum + 0.5;
    }
  }
}

int main(int argc, char* argv[])
{
  cv::Mat input;
  cv::Mat output;
  cv::Mat output_ref;
  cv::Mat output_cl;
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

  int rows, cols, steps;
  rows = input.rows;
  cols = input.cols;
  steps = input.step;
  unsigned char src[rows * cols], dst[rows * cols];

  int radius = (ksize.width - 1) / 2;

  output_cl = input.clone();
  output_ref = input.clone();
  output = input.clone();

  /* execute refrence codes */
  blur(input, output_ref, ksize, cv::Point(-1, -1), cv::BORDER_REFLECT_101);

  myblur(input, output, ksize);

  compareMat(output_ref, output);

#if 1
  int inputSize = input.rows * input.step;
  float fweight;

  // Load the kernel source code into the array source_str
  FILE *fp;
  char *source_str;
  size_t source_size;

  fp = fopen("blur.cl", "r");
  if (!fp) {
      fprintf(stderr, "Failed to load kernel.\n");
      exit(1);
  }
  source_str = (char*)malloc(MAX_SOURCE_SIZE);
  source_size = fread( source_str, 1, MAX_SOURCE_SIZE, fp);
  fclose( fp );

  // Get platform and device information
  cl_platform_id platform_id = NULL;
  cl_device_id device_id = NULL;   
  cl_uint ret_num_devices;
  cl_uint ret_num_platforms;
  cl_int ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
  ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_DEFAULT, 1, 
          &device_id, &ret_num_devices);

  // Create an OpenCL context
  cl_context context = clCreateContext( NULL, 1, &device_id, NULL, NULL, &ret);

  // Create a command queue
  cl_command_queue command_queue = clCreateCommandQueue(context, device_id,
      CL_QUEUE_PROFILING_ENABLE, &ret);

  // Create memory buffers on the device for each vector 
  cl_mem inputDev = clCreateBuffer(context, CL_MEM_READ_ONLY, 
          inputSize * sizeof(unsigned char), NULL, &ret);
  cl_mem dstDev = clCreateBuffer(context, CL_MEM_READ_WRITE, 
          inputSize * sizeof(unsigned char), NULL, &ret);

  // Copy the lists A and B to their respective memory buffers
  ret = clEnqueueWriteBuffer(command_queue, inputDev, CL_TRUE, 0,
          inputSize * sizeof(unsigned char), input.data, 0, NULL, NULL);


  // Create a program from the kernel source
  cl_program program = clCreateProgramWithSource(context, 1, 
          (const char **)&source_str, (const size_t *)&source_size, &ret);

  // Build the program
  ret = clBuildProgram(program, 1, &device_id, "-g -s \"blur_cl\"", NULL, NULL);

  if (ret != CL_SUCCESS) {
    size_t len;
    char buffer[2048];
    printf("Error: Failed to build program executable!\n");
    clGetProgramBuildInfo(program, device_id, 
    CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
    printf("%s\n", buffer);
    exit(1);
  }

  // Create the OpenCL kernel
  cl_kernel kernel = clCreateKernel(program, "blur_cl_local", &ret);

  if(ret != CL_SUCCESS) {
    cout << "Error: clCreateKernel" << endl;
    return 0;
  }


  fweight = (float)1 / (ksize.width * ksize.height);
  ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&inputDev);
  ret |= clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&dstDev);
  ret |= clSetKernelArg(kernel, 2, sizeof(int), (void *)&rows);
  ret |= clSetKernelArg(kernel, 3, sizeof(int), (void *)&cols);
  ret |= clSetKernelArg(kernel, 4, sizeof(int), (void *)&steps);

  /* be careful! it's unworked! */
  //ret |= clSetKernelArg(kernel, 4, sizeof(int), (void *)&(input.step));
  ret |= clSetKernelArg(kernel, 5, sizeof(int), (void *)&radius);
  ret |= clSetKernelArg(kernel, 6, sizeof(int), (void *)&fweight);
  ret |= clSetKernelArg(kernel, 7, 10000 * sizeof(unsigned char), NULL);

  if(ret != CL_SUCCESS) {
    cout << "Error: clSetKernelArg" << endl;
    return 0;
  }

  // Execute the OpenCL kernel on the list
  size_t global_item_size_2d[2];
  global_item_size_2d[0] = input.rows / 16 * 16;
  global_item_size_2d[1] = input.step / 16 * 16;
  cout << global_item_size_2d[0]<< ":" << global_item_size_2d[1] << endl;
  size_t local_item_size_2d[2];
  local_item_size_2d[0] = 16;
  local_item_size_2d[1] = 16;

  ret = clEnqueueNDRangeKernel(command_queue, kernel, 2, NULL, 
         global_item_size_2d, local_item_size_2d, 0, NULL, NULL);

  cout << "clEnqueueNDRangeKernel:" << ret << endl;

  for(int i = 0; i < inputSize; i++) {
    output_cl.data[i] = 0;
  }

  // Read the memory buffer C on the device to the local variable C
  ret = clEnqueueReadBuffer(command_queue, dstDev, CL_TRUE, 0, 
          inputSize * sizeof(unsigned char), output_cl.data, 0, NULL, NULL);

  cv::imshow("", output_cl);
  cv::waitKey(0);

  // Clean up
  ret = clFlush(command_queue);
  ret = clFinish(command_queue);
  ret = clReleaseKernel(kernel);
  ret = clReleaseProgram(program);
  ret = clReleaseMemObject(inputDev);
  ret = clReleaseMemObject(dstDev);
  ret = clReleaseCommandQueue(command_queue);
  ret = clReleaseContext(context);

  compareMat(output_ref, output_cl);
  //compareMat(input, output_cl);

  return 0;
#endif
}
