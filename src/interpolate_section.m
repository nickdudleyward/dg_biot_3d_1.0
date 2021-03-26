function vv = interpolate_section(h5_HW, h5_field, fld, x_samp, y_samp, z_samp)

% read dual mesh data and interpolation coefficients for the centroid of
% an element saved during mesh setup

centroid = h5read(h5_HW, '/centroid');
dual_tri = double(h5read(h5_HW, '/dual_tri'));
cen_interp = h5read(h5_HW, '/cen_interp');

% read field data for interpolation. Saved as a matrix where each column
% (in Fortran conventions) represents one field; select column specified in
% arguments and reshape it so it's indexed by local node numbers (row) and
% global element number (column)

field = h5read(h5_field, '/field');
field = reshape(field(:,fld), numel(cen_interp), []);

% now we can interpolate the field at the element centroids by multiplying 
% on the left by cen_interp (which is a row vector)

cen_field = cen_interp*reshape(field, NP, K);

% form all (x, y, z) triples from x_samp, y_samp, z_samp and save them 
% as rows of the matrix vv

[xx, yy, zz] = meshgrid(x_samp, y_samp, z_samp);
vv = [xx(:), yy(:), zz(:)];

% find tetrahedron numbers and barycentric coordinates of the elements
% of vv

[t, bary] = tsearchn(centroid, dual_tri, vv);

% identify sample points not in convex hull of dual mesh (i.e. t holds NaN)
% and mask out of t, bary and vv. Add an extra column to vv to hold 
% interpolated values of field

mask = ~isnan(t);
num_samples = nnz(mask);
t = t(mask);
bary = bary(mask, :);
vv = [ vv(mask, :), zeros(num_samples, 1) ];

% interpolate field at sample points, storing in the empty fourth column
% of vv created above
% Here, for the nth point, call it v:
% t(n) is the element of the dual mesh containing v, call it e
% tri(t(n),:) are the node numbers of the vertices of e
% field(tri(t(n),:)) are the field values at the vertices of e
% bary(n,:) are the barycentric coordinates of v wrt e
% the dot product of bary(...) and field(...) is the linear interpolation
% of field at v

for n=1:numel(t)
    vv(n, 4) = bary(n,:) * (cen_field(dual_tri(t(n),:)))';
end
