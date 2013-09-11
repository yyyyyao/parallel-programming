#ifndef __COMMON_H__
#define __COMMON_H__
#include <stdio.h>

#define dassert(err) do { \
  if(err < 0) { \
    printf("err:%d Error Occured@%s:%d:%s\n", err, __FILE__, __LINE__, __func__); \
    exit(0); \
  } \
} while(0)

#endif
