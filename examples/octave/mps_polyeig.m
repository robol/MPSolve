%% -*- texinfo -*-
%% @deftypefn {Function File} {@var{LAMBDA} =} mps_polyeig (@var{P0}, @var{P1}, ..., @var{Pn}, @var{alg} = 's')
%% @cindex root finding of a matrix polynomial
%% Compute the eigenvalues of the matrix polynomial p(z) given by
%% @tex
%% $$
%%  p(z) = P_0 + zP_1 + \\ldots + + z^{n}P_n 
%% $$ 
%% @end tex
%% @ifnottex
%% @example
%%         p(z) = P_0 + zP_1 + ... + z^nP_n 
%% @end example
%% @end ifnottex
%% and return a vector with the generalized eigenvalues.
%% The optional variable @var{alg} can be set to \"s\" or \"u\" to select 
%% the secular algorithm or the standard MPSolve algorithm. The default value 
%% is "s"
%% @end deftypefn
function LAMBDA = mps_polyeig(varargin)
  degree = length (varargin) - 1; 

  % In case the last parameter is a string treat it as the optional
  % algorithm parameter and lower the degree. 
  if (ischar (varargin{degree}))
    degree = degree - 1; 
    P = varargin; 
  else
    P = cell(1, degree + 2); 
    for i = 1:degree+1
      P{i} = varargin{i}; 
    endfor
    P{degree+2} = 's'; 
  endif

  % Warn if there are not a sufficient number of coefficients to
  % proceed. 
  if (degree <= 0)
    error ("mps_polyeig: please specify at least two coefficients for mps_polyeig! "); 
  endif

  % Handle a lot of special cases: 
  % 1) If the degree is 1 and we have that the linear term is well
  % conditioned transform the problem in a standard eigenvalue
  % problem. 
  if (degree == 1 && cond (varargin{1}) <= 1e4 * max(size(P{2})))
    P{1} = - P{2} \ P{1}; 
    P{2} = - eye (size (P{1})); 

    % Special code meaning that the problem is already Hessenberg. 
    P{degree+2} = 'h'; 

    % Take the problem in Hessenberg form. 
    [Q,H] = hess(P{1});
    LAMBDA = mps_polyeig_impl (H, - P{2}, [ P{3} 'h' ]);

  else
    LAMBDA = mps_polyeig_impl (varargin{:});
  endif

endfunction
