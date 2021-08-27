//
// File: fce_types.h
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 27-Feb-2017 22:42:17
//
#ifndef FCE_TYPES_H
#define FCE_TYPES_H

// Include Files
#include "rtwtypes.h"

// Type Definitions
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 //struct_emxArray_real_T

#ifndef struct_emxArray_uint8_T
#define struct_emxArray_uint8_T

struct emxArray_uint8_T
{
  unsigned char *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 //struct_emxArray_uint8_T
#endif

//
// File trailer for fce_types.h
//
// [EOF]
//
