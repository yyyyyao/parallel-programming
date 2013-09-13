/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2013年08月28日 17時26分14秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <iostream>
#include <math.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <CL/cl.h>

#define MAX_SOURCE_SIZE (0x100000)
//#define NO_SUBSAMPLING 1
#define INPUT_BMP 1

using namespace cv;
using namespace std;

int main(int argc, char const* argv[])
{
  int imgHeight = atoi(argv[1]);
  int imgWidth = atoi(argv[2]);
  int imgChannel = atoi(argv[3]);
  const char* fileName = argv[4];

  int imgSize = imgHeight * imgWidth;
  int imgMemSize = imgSize * 3;
#ifdef NO_SUBSAMPLING
  int halfImgHeight = imgHeight;
  int halfImgWidth = imgWidth;
#else
  int halfImgHeight = ceil(imgHeight / 2);
  int halfImgWidth = ceil(imgWidth / 2);
#endif

  unsigned char* pixel;
  pixel = (unsigned char*)malloc(sizeof(unsigned char) * imgMemSize);

  FILE* f = fopen(fileName, "rb");
  if(!f) {
    printf("file can't read!\n");
    return -1;
  }
  fread(pixel, sizeof(unsigned char), imgMemSize, f);
  fclose(f);

  unsigned char comp[imgChannel][imgSize];

  unsigned char*p = pixel;
  for(int i = 0; i < imgHeight * imgWidth; i++) {
    /* move to range 16-235 */
#ifdef INPUT_BMP
    comp[2][i] = (unsigned char)*(p++) * 0.8588 + 16;
    comp[1][i] = (unsigned char)*(p++) * 0.8588 + 16;
    comp[0][i] = (unsigned char)*(p++) * 0.8588 + 16;
#else
    comp[0][i] = (unsigned char)*(p++) * 0.8588 + 16;
    comp[1][i] = (unsigned char)*(p++) * 0.8588 + 16;
    comp[2][i] = (unsigned char)*(p++) * 0.8588 + 16;
#endif
  }

  cl_platform_id* platformId = NULL;
  cl_device_id deviceId = NULL;

  cl_uint retNumDevices;
  cl_uint retNumPlatforms;
  cl_int ret;

  ret = clGetPlatformIDs(1, NULL, &retNumPlatforms);
  if(ret != CL_SUCCESS) {
    printf("clGetPlatformIDS error!\n");
    return -1;
  }

  platformId = (cl_platform_id*)malloc(
      sizeof(cl_platform_id) * (unsigned int)retNumPlatforms);

  ret = clGetPlatformIDs(retNumPlatforms, platformId, &retNumPlatforms);
  if(ret != CL_SUCCESS) {
    printf("clGetPlatformIDS error!\n");
    return -1;
  }

  ret = clGetDeviceIDs(platformId[0], CL_DEVICE_TYPE_CPU, 1,
      &deviceId, &retNumDevices);
  if(ret != CL_SUCCESS) {
    printf("clGetDeviceIDs error!\n");
    return -1;
  }

  cl_context context = clCreateContext( NULL, 1, &deviceId, NULL, NULL, &ret);
  if(ret != CL_SUCCESS) {
    printf("clCreateContext error!\n");
    return -1;
  }

  cl_command_queue commandQueue = clCreateCommandQueue(context, deviceId,
      CL_QUEUE_PROFILING_ENABLE, &ret);
  if(ret != CL_SUCCESS) {
    printf("clCreateCommandQueue error!\n");
    return -1;
  }

  cl_mem srcMem = clCreateBuffer(context, CL_MEM_READ_WRITE,
      sizeof(unsigned char) * imgHeight * imgWidth * 3, NULL, &ret);
  if(ret != CL_SUCCESS) {
    printf("clCreateBuffer error!\n");
    return -1;
  }

  cl_mem dstMem = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
      sizeof(unsigned char) * imgHeight * imgWidth * 3, NULL, &ret);
  if(ret != CL_SUCCESS) {
    printf("clCreateBuffer error!\n");
    return -1;
  }
  cl_mem dstYMem = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
      sizeof(unsigned char) * imgHeight * imgWidth, NULL, &ret);
  if(ret != CL_SUCCESS) {
    printf("clCreateBuffer error!\n");
    return -1;
  }

  cl_mem dstCbMem = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
      sizeof(unsigned char) * halfImgHeight * halfImgWidth, NULL, &ret);
  if(ret != CL_SUCCESS) {
    printf("clCreateBuffer error!\n");
    return -1;
  }
  cl_mem dstCrMem = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
      sizeof(unsigned char) * halfImgHeight * halfImgWidth, NULL, &ret);
  if(ret != CL_SUCCESS) {
    printf("clCreateBuffer error!\n");
    return -1;
  }

  ret = clEnqueueWriteBuffer(commandQueue, srcMem, CL_TRUE, 0,
      sizeof(unsigned char) * imgHeight * imgWidth * 3, comp, 0, NULL, NULL);
  if(ret != CL_SUCCESS) {
    printf("clEnqueueWriteBuffer error!\n");
    return -1;
  }

  FILE* fp;
  char* sourceStr;
  size_t sourceSize;

  fp = fopen("myKernel.cl", "r");
  if(!fp) {
    printf("Failed to load kernel\n");
    return -1;
  }
  sourceStr = (char*)malloc(MAX_SOURCE_SIZE);
  sourceSize = fread(sourceStr, 1, MAX_SOURCE_SIZE, fp);
  fclose(fp);

  cl_program program = clCreateProgramWithSource(context, 1, (const char**)&sourceStr,
      (const size_t*)&sourceSize, &ret);
  if(program == NULL || ret != CL_SUCCESS) {
    printf("clCreateProgramWithSource error!\n");
    return -1;
  }

  ret = clBuildProgram(program, 1, &deviceId,
      "-g -s \"/home/yao/works/jpeg-origin/myKernel.cl\"", NULL, NULL);
  if(ret != CL_SUCCESS) {
    printf("clBuildProgram error!\n");
    char buffer[10240];
    clGetProgramBuildInfo(program, deviceId, CL_PROGRAM_BUILD_LOG,
        sizeof(buffer), buffer, NULL);
    printf("compile error!:\n%s\n", buffer);
    return -1;
  }

  cl_kernel kernel = clCreateKernel(program, "colConv_function", &ret);
  if(kernel == NULL || ret != CL_SUCCESS) {
    printf("clCreateKernel error!\n");
    return -1;
  }

  clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&srcMem);
  clSetKernelArg(kernel, 1, sizeof(cl_mem), (void*)&dstMem);
  clSetKernelArg(kernel, 2, sizeof(int), (void*)&imgWidth);
  clSetKernelArg(kernel, 3, sizeof(int), (void*)&imgHeight);
  clSetKernelArg(kernel, 4, sizeof(cl_mem), (void*)&dstYMem);
  clSetKernelArg(kernel, 5, sizeof(cl_mem), (void*)&dstCbMem);
  clSetKernelArg(kernel, 6, sizeof(cl_mem), (void*)&dstCrMem);

  size_t globalItemSize[2];
  globalItemSize[0] = 256;
  globalItemSize[1] = 256;
  size_t localItemSize[2];
  localItemSize[0] = 4;
  localItemSize[1] = 4;

  cl_event event;
  clEnqueueNDRangeKernel(commandQueue, kernel, 2, NULL, globalItemSize, localItemSize,
      0, NULL, &event);
  clWaitForEvents(1, &event);

  cl_ulong timeStart, timeEnd;
  double totalTime;
  ret = clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START,
      sizeof(cl_ulong), &timeStart, NULL);
  if(ret != CL_SUCCESS) {
    printf("clGetEventProfilingInfo error!\n");
    return -1;
  }
  ret = clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END,
      sizeof(cl_ulong), &timeEnd, NULL);
  if(ret != CL_SUCCESS) {
    printf("clGetEventProfilingInfo error!\n");
    return -1;
  }

  totalTime = timeEnd - timeStart;
  printf("Execution time[msec] %0.3f\n", totalTime / 1000000);

  unsigned char dstImageY[imgHeight * imgWidth];
  ret = clEnqueueReadBuffer(commandQueue, dstYMem, CL_TRUE, 0,
      imgHeight * imgWidth * sizeof(unsigned char), dstImageY, 0, NULL, NULL);
  if(ret != CL_SUCCESS) {
    printf("clEnqueueReadBuffer error!\n");
    return -1;
  }


  unsigned char dstImageCb[halfImgHeight * halfImgWidth];
  unsigned char dstImageCr[halfImgHeight * halfImgWidth];
  ret = clEnqueueReadBuffer(commandQueue, dstCbMem, CL_TRUE, 0,
      halfImgHeight * halfImgWidth * sizeof(unsigned char), dstImageCb, 0, NULL, NULL);
  if(ret != CL_SUCCESS) {
    printf("clEnqueueReadBuffer error!\n");
    return -1;
  }

  ret = clEnqueueReadBuffer(commandQueue, dstCrMem, CL_TRUE, 0,
      halfImgHeight * halfImgWidth * sizeof(unsigned char), dstImageCr, 0, NULL, NULL);
  if(ret != CL_SUCCESS) {
    printf("clEnqueueReadBuffer error!\n");
    return -1;
  }

  Mat dst(imgWidth, imgHeight, CV_8UC3);
  unsigned char* dstPtr = dst.ptr<unsigned char>();

#ifdef NO_SUBSAMPLING
  for(int i = 0; i < imgWidth; i++) {
    for(int j = 0; j < imgHeight; j++) {
      float y = dstImageY[i * imgWidth + j];
      float cb = dstImageCb[i * imgWidth + j];
      float cr = dstImageCr[i * imgWidth + j];
      dstPtr[i * imgWidth * 3 + j * 3] = (unsigned char) (y + 1.4 * (cr - 128));
      dstPtr[i * imgWidth * 3 + j * 3 + 1] = (unsigned char) (y - 0.344 * (cb - 128) - 0.714 * (cr - 128));
      dstPtr[i * imgWidth * 3 + j * 3 + 2] = (unsigned char) (y + 1.772 * (cb - 128));
    }
  }

#else
  for(int i = 0; i < imgWidth; i++) {
    for(int j = 0; j < imgHeight; j++) {
      float y = dstImageY[i * imgWidth + j];
      float cb = dstImageCb[(i / 2) * halfImgWidth + j / 2];
      float cr = dstImageCr[(i / 2) * halfImgWidth + j / 2];

      dstPtr[i * imgWidth * 3 + j * 3] = (unsigned char) (y + 1.4 * (cr - 128));
      dstPtr[i * imgWidth * 3 + j * 3 + 1] = (unsigned char) (y - 0.343 * (cb - 128) - 0.711 * (cr - 128));
      dstPtr[i * imgWidth * 3 + j * 3 + 2] = (unsigned char) (y + 1.765 * (cb - 128));
    }
  }
#endif

  imshow("dst", dst);
  waitKey(0);

  clFlush(commandQueue);
  clFinish(commandQueue);
  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseMemObject(srcMem);
  clReleaseMemObject(dstMem);
  clReleaseCommandQueue(commandQueue);
  clReleaseContext(context);

  return 0;

}
