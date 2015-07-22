% Generate a random example and try to solve it. 
n = 30; p = 3; a = randn(n,1); b = randn(n,1); c = randn(n,1);

x = rational_roots(a, b, c, p)

% Check if the residue is small. 
for i = 1 : length(x)
  fprintf ('Residue = %e\n', sum (c ./ (a + b * x(i)).^p));
end
