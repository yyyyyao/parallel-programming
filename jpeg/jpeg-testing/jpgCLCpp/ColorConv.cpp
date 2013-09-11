#include <stdlib.h>
#ifdef __APPLE__
#include<OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#include <string>
#include <stdio.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "include/jpgInternal.h"
#include "include/jpgCLInternal.h"
#include "ColorConv.h"

/* you need to fix this func. */
ColorConv::ColorConv(cv::Mat* inputMat, ClContext* cl) {
  m_cl = cl;
  m_image_x = inputMat->cols;
  m_image_y = inputMat->rows;
  m_image_half_x = m_image_x / 2;
  m_image_half_y = m_image_y / 2;
  m_inputPtr = (uint8*)inputMat->data;
  m_Y.create(m_image_x, m_image_y, CV_8UC1);
  m_YPtr = (uint8*)m_Y.data;
  m_Cb.create(m_image_half_x, m_image_half_y, CV_8UC1);
  m_CbPtr = (uint8*)m_Cb.data;
  m_Cr.create(m_image_half_x, m_image_half_y, CV_8UC1);
  m_CrPtr = (uint8*)m_Cr.data;
  m_sourceName = (char*)"kernel/ColorConv.cl";
  m_funcName = (char*)"ColorConv_SubSampling";
  m_absFilePath = (char*)"";
  m_kernelDimension = 2;
  m_globalItemSize[0] = m_image_x >> 1;
  m_globalItemSize[1] = m_image_y >> 1;
  m_localItemSize[0] = 2;
  m_localItemSize[1] = 2;

  buildProgram();
  setKernelArgs_Execute();
}

/* you need to fix this func. */
void ColorConv::setKernelArgs_Execute() {
  int ret;
  cl_event event;
  cl_ulong timeStart = 0, timeEnd = 0;

  /* Allocate Device Memory. */
  cl_mem srcMemObj = clCreateBuffer(m_cl->m_context, CL_MEM_READ_WRITE,
      sizeof(unsigned char) * m_image_x * m_image_y * 3, NULL, &ret);
  dassert(ret);
  cl_mem dstYMemObj = clCreateBuffer(m_cl->m_context, CL_MEM_WRITE_ONLY,
      sizeof(unsigned char) * m_image_x * m_image_y, NULL, &ret);
  dassert(ret);
  cl_mem dstCbMemObj = clCreateBuffer(m_cl->m_context, CL_MEM_WRITE_ONLY,
      sizeof(unsigned char) * m_image_half_x * m_image_half_y, NULL, &ret);
  dassert(ret);
  cl_mem dstCrMemObj = clCreateBuffer(m_cl->m_context, CL_MEM_WRITE_ONLY,
      sizeof(unsigned char) * m_image_half_x * m_image_half_y, NULL, &ret);
  dassert(ret);

  /* Write Data. */
  clEnqueueWriteBuffer(m_cl->m_commandQueue, srcMemObj, true, 0,
      sizeof(unsigned char) * m_image_x * m_image_y * 3, m_inputPtr, 0, NULL, NULL);
  dassert(ret);
  /* Set Kernel Args */
  ret = clSetKernelArg(m_kernel, 0, sizeof(cl_mem), (void *)&srcMemObj);
  dassert(ret);
  ret |= clSetKernelArg(m_kernel, 1, sizeof(int), (void *)&m_image_x);
  dassert(ret);
  ret |= clSetKernelArg(m_kernel, 2, sizeof(int), (void *)&m_image_y);
  dassert(ret);
  ret |= clSetKernelArg(m_kernel, 3, sizeof(cl_mem), (void *)&dstYMemObj);
  dassert(ret);
  ret |= clSetKernelArg(m_kernel, 4, sizeof(cl_mem), (void *)&dstCbMemObj);
  dassert(ret);
  ret |= clSetKernelArg(m_kernel, 5, sizeof(cl_mem), (void *)&dstCrMemObj);
  dassert(ret);

  ret = clEnqueueNDRangeKernel(m_cl->m_commandQueue, m_kernel, 2, NULL,
      m_globalItemSize, m_localItemSize, 0, NULL, &event);
  dassert(ret);

  clWaitForEvents(1, &event);
  clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START,
      sizeof(cl_ulong), &timeStart, NULL);
  clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END,
      sizeof(cl_ulong), &timeEnd, NULL);

  printf("Func:%s Execution time:%0.3f[msec]\n",m_funcName, (double)(timeEnd - timeStart) / 1000000);

  clEnqueueReadBuffer(m_cl->m_commandQueue, dstYMemObj, CL_TRUE, 0,
      m_image_x * m_image_y * sizeof(unsigned char), m_YPtr, 0, NULL, NULL);
  clEnqueueReadBuffer(m_cl->m_commandQueue, dstCbMemObj, CL_TRUE, 0,
      m_image_half_x * m_image_half_y * sizeof(unsigned char), m_CbPtr, 0, NULL, NULL);
  clEnqueueReadBuffer(m_cl->m_commandQueue, dstCrMemObj, CL_TRUE, 0,
      m_image_half_x * m_image_half_y * sizeof(unsigned char), m_CrPtr, 0, NULL, NULL);

  /* Release Kenel, Program, and Memory. */
  clReleaseKernel(m_kernel);
  clReleaseProgram(m_program);
  clReleaseMemObject(srcMemObj);
  clReleaseMemObject(dstYMemObj);
  clReleaseMemObject(dstCbMemObj);
  clReleaseMemObject(dstCrMemObj);

#if 1
  /* showImage */
  cv::namedWindow("Test", CV_WINDOW_AUTOSIZE);
  cv::imshow("Test", m_Cr);
  cv::waitKey(0);
#endif
}

/* you don't have to fix this func.
 * use m_sourceName/m_funcName */
void ColorConv::buildProgram() {
  FILE* fp;
  char* sourceStr;
  char clBuildOptions[BUF_SIZE];
  size_t sourceSize;
  int ret;

  fp = fopen(m_sourceName, "r");
  if(!fp) {
    fprintf(stderr, "Failed to load kernel.\n");
    exit(1);
  }

  sourceStr = (char*)malloc(sizeof(char) * MAX_SOURCE_SIZE);
  /* should be checked. */
  sourceSize = fread(sourceStr, 1, MAX_SOURCE_SIZE, fp);
  fclose(fp);


  m_program = clCreateProgramWithSource(m_cl->m_context, 1,
      (const char **)&sourceStr, (const size_t*)&sourceSize, &ret);
  dassert(ret);

  switch(m_cl->m_pType) {
    case Intel:
      sprintf(clBuildOptions, "-D INTEL_PLATFORM -g %s", m_absFilePath);
      break;
    case NVIDIA:
      sprintf(clBuildOptions, "-D NVIDIA_PLATFORM");
      break;
    default:
      sprintf(clBuildOptions, "");
      break;
  }

  ret = clBuildProgram(m_program, 1, &(m_cl->m_deviceId), clBuildOptions, NULL, NULL);
  if(ret != CL_SUCCESS) {
    char buffer[MAX_SOURCE_SIZE];
    clGetProgramBuildInfo(m_program, m_cl->m_deviceId, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, NULL);
    fprintf(stderr, "CL Compilation failed\n%s", buffer);
    exit(1);
  }

  m_kernel = clCreateKernel(m_program, m_funcName, &ret);
  dassert(ret);
}

void ColorConv::setGlobalWorkGroupSize(int dim, size_t* gWork) {
  if(dim < 1 || dim > 3) {
    printf("invalid value@dim:%d\n", dim);
    return;
  }
  for(int i = 0; i < dim; i++) {
    m_globalItemSize[i] = gWork[i];
  }
}

void ColorConv::setLocalWorkGroupSize(int dim, size_t* lWork) {
  if(dim < 1 || dim > 3) {
    printf("invalid value@dim:%d\n", dim);
    return;
  }
  for(int i = 0; i < dim; i++) {
    m_localItemSize[i] = lWork[i];
  }
}

ColorConv::~ColorConv() {

}
