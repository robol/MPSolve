%% MPS_ROOTS function that approximate the roots of a polynomial specified
%% by the vector v using MPSolve and the algorithm specified by alg. 
%% alg can be 's' or 'u', depending if the algorithm MPS_ALGORITHM_SECULAR_GA or
%% the algorithm MPS_ALGORITHM_SECULAR_UNISOLVE are desired. 
%%
%% Author: Leonardo Robol <leonardo.robol@sns.it>
%% Copyright: 2011-2013 Leonardo Robol <leonardo.robol@sns.it>
%% License: GPLv3 or higher
function x = mps_roots(v, alg)
    
  if min(size(v)) ~= 1 || strcmp(class(v(1)), 'string')
    error('The input must be a 1D vector');
  end

  if nargin <= 1
     alg = 's';
  end

  if isnumeric(v(1))
    x = mps_roots_double (v, alg);
  else
    is_vpa = strcmp (class(v(1)), 'sym');
    
    % If the input is given in VPAs, convert them to string.
    if is_vpa
      vv = cell(1,length (v));
      for i = 1 : length (v)
        vv{i} = char(v(i));
      end
      v = vv;
    end
    
    x = mps_roots_string (vv, alg);
  end
end
