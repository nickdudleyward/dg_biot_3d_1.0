function [elt, interp] = locate_element(g)

% From a matrix g representing a number of points (each row in a coordinate
% vector), find the containing element and the interpolation vector from
% the nodal values on that element to the field value at that point

Globals3D;

% find elements containing specified points and barycentric coordinates
% of these points within elements
% g is an (n x 3) matrix whose rows contain the points' coords

[elt, bary] = tsearchn([VX',VY',VZ'], EToV, g);

% elt is a column vector of n element numbers and bary is an
% (n x 4) matrix of barycentric coordinates

% Convert barycentrics to r, s, t coordinate system: H&W p.409
% rst is an (n x 3) matrix whose rows contains local element 
% coordinates of points in g

rst = bary*[[-1, -1, -1]; [1, -1, -1]; [-1, 1, -1]; [-1, -1, 1]];

% Vandermonde matrix for evaluating field at each point in g
% (from field modal values on each receiver element). V is (n x Np)
% for n points, i.e. each row holds the interpolation vector for a
% receiver
% Note that Vandermonde3D behaves very strangely if it gets rows instead
% of columns for its r, s, t arguments!
V0 = Vandermonde3D(N, rst(:,1), rst(:,2), rst(:,3));
% Interpolation matrix for evaluating field at each receiver point
% (from field nodal values on each receiver element): global invV maps 
% from nodal to modal values.
interp = V0*invV;
