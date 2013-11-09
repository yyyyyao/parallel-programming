#ifndef __CLCONTEXTHPP__
#define __CLCONTEXTHPP__

#include "common.h"
#include "include/jpgCLInternal.h"
#ifdef __APPLE__
#include<OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

class ClContext{
  public:
    /* for OpenCL variables */
    cl_platform_id* m_platformId;
    cl_device_id m_deviceId;
    cl_context m_context;
    cl_command_queue m_commandQueue;
    cl_uint m_platformNum;
    cl_uint m_deviceNum;
    platformType m_pType;

    ClContext(const char* pName);
    ClContext();
    ~ClContext();
};
#endif
