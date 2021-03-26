% Fix the axis ax ('x', 'y' or 'z') at value w. 
% On the first of the other two axis (in the order x, y, z) sample
% M points from u0 to u1.
% On the second of the other two axis (in the order x, y, z) sample
% N points from v0 to v1.
% Combine these with w to create an MxN mesh of 3D coordinates in the 
% plane ax = w
% Sample field number fld on these points

Globals3D;

ax = 'x'; w = 1.0; 

u_vals = linspace(u0, u1, M);
v_vals = linspace(v0, v1, N);

if strcmp(ax, 'x')
    x_mesh, y_mesh, zmesh = meshgrid(w, u_vals, v_vals);
elseif strcmp(ax, 'y')
    x_mesh, y_mesh, zmesh = meshgrid(u_vals, w, v_vals);
elseif  strcmp(ax, 'z')
    x_mesh, y_mesh, zmesh = meshgrid(u_vals, v_vals, w);
else
    error('Unrecognised axis %s\n', ax);
end

% The three _mesh are 3D arrays. Flatten them and assemble as an (M*N x 3)
% matrix where each row holds the coordinates of a mesh point.

xyz_mesh = [x_mesh(:), y_mesh(:), z_mesh(:)];

% For each mesh point, find the element containing it and the barycentric
% coordinates of the mesh point within the containing element

[elt_mesh, bary_mesh] = tsearchn([VX',VY',VZ'], EToV, xyz_mesh);

% elt is a column vector of M*N element numbers and bary is an
% (M*N x 4) matrix of barycentric coordinates

% Convert barycentrics to r, s, t coordinate system: H&W p.409
% rst_rec is a (3 x n) matrix whose columns contains local element 
% coordinates of receivers

rst_mesh = [-1,1,-1,-1; -1,-1,1,-1;-1,-1,-1,1]*bary';

% Vandermonde matrix for evaluating field at each receiver point
% (from field modal values on each receiver element). V_rec is (M*N x Np)

V_mesh = Vandermonde3D(N, rst_mesh(1,:), rst_mesh(2,:), rst_mesh(3,:));

% Interpolation matrix for evaluating field at each mesh point
% (from field nodal values on each receiver element): global Vinv maps 
% from nodal to modal values.

interp_mesh = V_mesh*invV;

% interpolate field fld at each point in the mesh

w_vals = zeros(M*N);
for j=1:M*N
    w_vals(j) = interp_mesh(j, :)*field(:, elt_mesh(j), fld);
end

% create a 2D mesh on the u, v axes for plotting and reshape the values
% in w_vals consistent with this

[u_mesh, v_mesh] = meshgrid(u_vals, v_vals);
w_vals = reshape(w_vals, M, N);












% Slice through a mesh

Nx=1; Ny=2; Nz=3;

% Polar coordinate conventions: for a point (x,y,z) in Cartesian 
% coordinates and (r, theta, phi) in spherical mpolar cordinates:
% theta from -pi to pi is the azimuth angle, i.e. the angle from the 
% positive x axis of the projection of the vector into the (x,y) plane
% phi from 0 to pi is the zenith angle, i.e. the angle from the positive 
% z axis of the vector
% We have
% x = r*cos(theta)*sin(phi)
% y = r*cos(theta)*sin(phi)
% z = r*cos(phi)
% and
% r     = sqrt(x^2+y^2+z^2)
% theta = atan2(y,x)
% phi   = acos(z/r)

% We're given a normal vector N=(Nx,Ny,Nz)
% Start off by converting to spherical polar coordinates

Nr = sqrt(Nx^2+Ny^2+Nz^2);
Ntheta = atan2(Ny,Nx);
Nphi = acos(Nz/Nr);

% Now find an orthonormal basis containing the unit 
% vector[1, Ntheta, Npi]. These are normalised copies of 
% d(r,theta,phi)/dtheta and d(r,theta,phi)/dphi

b = [Nx; Ny; Nz]/Nr;
btheta = [-sin(Ntheta); cos(Ntheta); 0];
bphi = [cos(Ntheta)*cos(Nphi); sin(Ntheta)*cos(Nphi); -sin(Nphi)];

% Assemble these as colummns of a matrix U.

U = [bphi, btheta, b];

% U is unitary, orientation preserving and maps [0,0,1] to the unit vector
% (Nx,Ny,Nz)/Nr
% Check: all these should be zero
U*U'-eye(3)
det(U)-1
U*[0;0;1]-[Nx;Ny;Nz]/Nr

% Construct a mesh in the [u,v,0] plane and map it to [x,y,z] space using U
