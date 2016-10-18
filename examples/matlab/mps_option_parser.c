#include "mps_option_parser.h"
#include <string.h>

_mps_matlab_options
mps_parse_matlab_options (const mxArray * optionStruct)
{
  _mps_matlab_options options = {
    MPS_ALGORITHM_SECULAR_GA,
    MPS_OUTPUT_GOAL_APPROXIMATE,
    16,
    false
  };

  if (!optionStruct)
    return options;

  if (mxIsChar (optionStruct))
    {
      mxChar * value = mxGetChars (optionStruct);
      if (strcmp ((char*) value, "s") == 0)
	options.algorithm = MPS_ALGORITHM_SECULAR_GA;
      else if (strcmp ((char*) value, "u") == 0)
	options.algorithm = MPS_ALGORITHM_STANDARD_MPSOLVE;
      else
	mexErrMsgTxt ("Invalid value specified for the algorithm. Only 'u' or 's' are allowed.\n");      
    }
  else if (! mxIsStruct (optionStruct))
    mexErrMsgTxt ("Only chars and struct values are allowed as MPSolve options\n");

  if (optionStruct)
    {
      int nFields = mxGetNumberOfFields (optionStruct);
      int i;

      for (i = 0; i < nFields; i++)
	{
	  const char * optionName = mxGetFieldNameByNumber (optionStruct, i);
	  mxArray * field = mxGetFieldByNumber (optionStruct, 0, i);

	  if (strcmp (optionName, "algorithm") == 0)
	    {
	      if (! mxIsChar (field))
		mexErrMsgTxt ("Please specify only 'u' or 's' for the algorithm property\n");
	      else
		{	      
		  mxChar * value = mxGetChars (field);
		  if (strcmp ((char*) value, "s") == 0)
		    options.algorithm = MPS_ALGORITHM_SECULAR_GA;
		  else if (strcmp ((char*) value, "u") == 0)
		    options.algorithm = MPS_ALGORITHM_STANDARD_MPSOLVE;
		  else
		    mexErrMsgTxt ("Invalid value specified for the property: 'algorithm'. Only 'u' or 's' are allowed.\n");
		}
	    }
	  else if (strcmp (optionName, "digits") == 0)
	    {
	      if (! mxIsNumeric (field))
		mexErrMsgTxt ("Please specify a positive integer for the digits property\n");
	      else
		{
		  double digits = mxGetScalar (field);
	      
		  if (digits <= 0)
		    mexErrMsgTxt ("Please specify a positive integer for the digits property\n");
		  else
		    options.digits = (int) digits;
		}
	    }
	  else if (strcmp (optionName, "goal") == 0)
	    {
	      if (! mxIsChar (field))
		mexErrMsgTxt ("Please specify only 'a' or 'i' as value for the goal property\n");
	      else
		{
		  mxChar * value = mxGetChars (field);
		  if (strcmp ((char*) value, "i") == 0)
		    options.goal = MPS_OUTPUT_GOAL_ISOLATE;
		  else if (strcmp ((char*) value, "a") == 0)
		    options.goal = MPS_OUTPUT_GOAL_APPROXIMATE;
		  else
		    mexErrMsgTxt ("Please specify only 'a' or 'i' as value for the goal property\n");		    
		}
	    }
	  else if (strcmp (optionName, "radius") == 0)
	    {
	      if (! mxIsLogicalScalar (field))
		mexErrMsgTxt ("Please specify 'true' or 'false' as a value for the radius property\n");
	      else
		{
		  options.radius = *mxGetLogicals (field);
		}
	    }
	  else
	    {
	      char * buffer = (char*) mxMalloc (sizeof (char) * (strlen ("The property '' is invalid\n") + strlen ((char*) optionName) + 1));
	      sprintf (buffer, "The property '%s' is invalid\n", optionName);
	      mexErrMsgTxt (buffer);
	      mxFree (buffer);
	    }
	}
    }

  return options;
}


