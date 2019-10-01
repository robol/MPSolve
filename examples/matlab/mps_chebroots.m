function [x, r] = mps_chebroots(v)
%%MPS_CHEBROOTS Compute roots of polynomials in the Chebyshev basis. 
%
% X = MPS_CHEBROOTS(V) computes the roots of the polynomial 
%
%  V(1) * T_N(X) + ... + V(end) * T_0(X), 
%
% where T_j(X) is the degree j Chebyshev polynomial of the first kind. 

opts = struct;

opts.chebyshev = true;
opts.radius = nargout > 1;

if opts.radius
	[x, r] = mps_roots_double(v, opts);
else
	x = mps_roots_double(v, opts);
end
 
end
