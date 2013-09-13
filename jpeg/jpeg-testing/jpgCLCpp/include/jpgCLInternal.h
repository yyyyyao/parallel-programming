#ifndef __JPGCLINTERNALH__
#define __JPGCLINTERNALH__
//typedef enum {NVIDIA, Intel, AMD, OTHERS} platformType;
typedef enum {Intel, NVIDIA, AMD, OTHERS} platformType;
//const platformType platformList[] = {NVIDIA, Intel, AMD, OTHERS};
const platformType platformList[] = {Intel, NVIDIA, AMD, OTHERS};
#define MAX_SOURCE_SIZE 16384
#define BUF_SIZE 2048
#endif

