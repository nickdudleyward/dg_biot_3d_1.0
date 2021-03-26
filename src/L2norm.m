function n = L2norm(u)

% L^2 norm across the whole mesh of the (element) node values u

Globals3D

% invV comverts the nodal values into modal, i.e. coefficients wrt
% an orthonormal basis on the reference element. Multiplication by J 
% corrects for the measure of the actual element.

n = sqrt(sum(sum(((invV*u).^2).*J)));
