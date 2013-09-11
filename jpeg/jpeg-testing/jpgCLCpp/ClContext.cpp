#include <iostream>
#include <typeinfo>
#include <string>
#include <stdio.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "include/jpgInternal.h"
#include "include/jpgCLInternal.h"
#include "include/ClContext.h"

#define PLATFORM_SIZE 10

#ifdef __APPLE__
#include<OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

ClContext::ClContext(const char* pName) {
  int ret;
  clGetPlatformIDs(1, NULL, &m_platformNum);
  m_platformId = (cl_platform_id*)malloc(sizeof(cl_platform_id) * m_platformNum);

  clGetPlatformIDs(m_platformNum, m_platformId, &m_platformNum);

  /* choose platform. Intel or NVIDIA */
  char platformName[PLATFORM_SIZE][BUF_SIZE];

  /* show & get platform names.
   * this code is incomplete, because
   * it depends on order of clGetPlatformInfo*/
  for(int i = 0; i < (int)m_platformNum; i++) {
    ret = clGetPlatformInfo(m_platformId[i], CL_PLATFORM_VENDOR,
        BUF_SIZE, &platformName, NULL);
    dassert(ret);
    printf("%d: %s\n", i, platformName[i]);
  }
  for(int i = 0; i < (int)m_platformNum; i++) {
    if(strstr(platformName[i], pName) != NULL) {
      printf("set platform:%s\n", platformName[i]);
      m_pType = platformList[i];
      ret = clGetDeviceIDs(m_platformId[i], CL_DEVICE_TYPE_DEFAULT, 1,
          &m_deviceId, &m_deviceNum); 
      dassert(ret);
      break;
    }
  }

  m_context = clCreateContext(NULL, 1, &m_deviceId, NULL, NULL, &ret);
  dassert(ret);

  m_commandQueue = clCreateCommandQueue(m_context, m_deviceId,
      CL_QUEUE_PROFILING_ENABLE, &ret);
  dassert(ret);
}

#if 0
ClContext::ClContext() {
  printf("hogehoge\n");
}
#endif
ClContext::~ClContext() {
  free(m_platformId);
  clFlush(m_commandQueue);
  clFinish(m_commandQueue);
  clReleaseCommandQueue(m_commandQueue);
  clReleaseContext(m_context);
}
