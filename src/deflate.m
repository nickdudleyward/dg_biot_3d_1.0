% A should be a vector of node numbers of local size n, i.e. a
% concatenation of blocks of length n of the form 
% (k-1)*n+1, (k-1)*n+2, ... k*n
% Each such block is replaced by the value k
% Intended use: from a list of e.g. face nodes, find th corresponding face
% numbers.
function B = deflate(A, n)

B = A(1:n:end);
B = ((B-1)/n)+1;

