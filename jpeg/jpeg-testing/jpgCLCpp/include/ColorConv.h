#ifndef __COLORCONV_H__
#define __COLORCONV_H__

#include "include/common.h"
#include "include/ClContext.h"
class ColorConv {
  public:
    ColorConv(cv::Mat*, ClContext*); 
    ~ColorConv();
    void execute();

    void setGlobalWorkGroupSize(int , size_t*);
    void setLocalWorkGroupSize(int , size_t*);
  private:
    void buildProgram();
    void setKernelArgs_Execute();
    /* private datas. */
    int m_image_y;
    int m_image_x;
    int m_image_half_x;
    int m_image_half_y;
    unsigned char* m_inputPtr;
    cv::Mat m_Y;
    cv::Mat m_Cb;
    cv::Mat m_Cr;
    unsigned char* m_YPtr;
    unsigned char* m_CbPtr;
    unsigned char* m_CrPtr;

   /* for OpenCL variables */
   ClContext* m_cl; 
   cl_program m_program;
   cl_kernel m_kernel;
   cl_mem m_srcMem;
   cl_mem m_dstYMem;
   cl_mem m_dstCbMem;
   cl_mem m_dstCrMem;
   char* m_sourceName;
   char* m_funcName;
   char* m_absFilePath; /* for Intel OpenCL Debugger. */
   unsigned int m_kernelDimension;
   size_t m_globalItemSize[2]; /* max dimension */
   size_t m_localItemSize[2];
};
#endif
