#include <mex.h>
#include <mps/mps.h>

typedef struct {
  mps_algorithm algorithm;
  mps_output_goal goal;
  int digits;
} _mps_matlab_options;

_mps_matlab_options mps_parse_matlab_options (const mxArray * optionStruct);


    
  
