function plot_section_3d(xval, yval, zval, f, num_points, varargin)

% Plots a section through a function defined on the points of a 3D mesh.
%
% Derived from Python code in FEM.plot_section_3d but uses a DG-style
% mesh and higher-order interpolation, for which it requires a
% Hesthaven-Warburton 3D mesh structure.
% 
% Exactly one of xval, yval, zval should be a scalar; the other two should
% each be a vector of length 2. The scalar indicates an x, y, or z
% coordinate to fix, defining a plane in 3-dimensional space. The two
% vectors of length 2 indicate minimal and maximal values of x, y, z
% coordinates, defining a rectangle in that plane. f should be a vector of
% values on the nodes of the Hesthave-Warburton style mesh, the same
% shape as e.g. the global x, y, z matrices. num_points  should be a 
% positive integer.
%
% The data in f is interpolated across an num_points by num_points array of 
% equally-spaced samples in the rectangle defined by xval, yval, zval. 
% It is then plotted in one of three ways, according to a string passed
% as an options extra argument
% 
% 'heat' (default) plots a rectangular array of cells, coloured according 
% to the value of the function in the centre of each cell
% 'contour' plots contours. HOW TO SPECIFY COUNTOUR LEVELS?
% 'gradient' plots an arrow at each point, indicating the gradient vector 
% at that point.
% Extra options for 'gradient' (ignored in other styles)
% 'unit_gradient' gives each vector the same length
% 'reverse_gradient' changes the direction of each vector

%% STEP 1: process arguments and find mesh onto which field data will ne
%% interpolated

Globals3D;

% Boolean options

bools = { 'unit_gradient', 'reverse_gradient' };

% Style options

styles = { 'heat', 'gradient', 'contour' };

% Initialise all options to false

for s = [bools, styles]
    opt.(s{1}) = false;
end

for j=1:numel(varargin)
    arg = varargin{j};
    if ~isfield(opt, arg)
        error('Unrecognised option %s', arg)
    end
    opt.(arg) = true;
end

% count number of style options

num_styles = 0;
for s = styles
    if opt.(s{1})
        num_styles = num_styles + 1;
    end
end

% if there are no style options, use 'heat' as the default
% if there are more than one style option, trigger an error

if num_styles == 0
    opt.heat = true;
elseif num_styles > 1
    error('More than one plotting style specified');
end

% Looks daft, but it leaves the opportunity for expansion!
% For heat map, we sample in the middle of each rectangle;
% for other plots, we sample on the corners

%centre_sampling = opt.heat;
centre_sampling = false;

% flatten xval, yval, zval in case they're an odd shape

xval = xval(:); yval = yval(:); zval = zval(:); 
    
% each of xval, yval, zval length 1 (a point) or 2 (a range)

if numel(xval) < 1 || numel(xval) > 2 || numel(yval) < 1 || numel(yval) > 2 || numel(zval) < 1 || numel(zval) > 2
    error('x, y, z values must have length 1 (a point) or 2 (a range)');
end

% for a section, need one point and two ranges

if numel(xval) + numel(yval) + numel(zval) == 4
    sample_dim = 1;
elseif numel(xval) + numel(yval) + numel(zval) == 5
    sample_dim = 2;
else
    error('x, y, z values invalid');
end

% convert the ranges into arrays of N sample points. If centre_sampling
% is True, these are the centres of subintervals of equal width
% covering the sampling interval; if False, they are at endpoints

x_offset = 0.0; y_offset = 0.0; z_offset = 0.0;

if numel(xval) == 2
    if centre_sampling
        x_offset = 0.5*(valx(2)-xval(1))/(num_points+1);
    end
    xx = linspace(xval(1)+x_offset, xval(2)-x_offset, num_points);
end
if numel(yval) == 2
    if centre_sampling
        y_offset = 0.5*(yval(2)-yval(1))/(num_points+1);
    end
    yy = linspace(yval(1)+y_offset, yval(2)-y_offset, num_points);
end
if numel(zval) == 2
    if centre_sampling
        z_offset = 0.5*(zval(2)-zval(1))/(num_points+1);
    end
    zz = linspace(zval(1)+z_offset, zval(2)-z_offset, num_points);
end

% Set up arrays X, Y, Z of coordinates of sample points.
%
% For a 1-dimensional plot, one of these is a linspace between the
% two sample points, the other two are constant of the same size
%
% For a two-dimensional plot, meshgrid the two lispaces from the ranges 
% together to give pairs of sample points and makethe other constant of
% the same size. Record the maximum and minimum values of the two ranges 
% as umin, umax (horizontal axis) and vmin, vmax (vertical axis)

if sample_dim == 1
    if numel(xval) == 2
        X = xx;
        Y = yval*ones(size(X));
        Z = zval*ones(size(X));
        UU = X;
        VV = [];
    end
    if numel(yval) == 2
        Y = yy;
        X = xval*ones(size(Y));
        Z = zval*ones(size(Y));
        UU = Y;
        VV = [];
    end
    if numel(zval) == 2
        Z = zz;
        X = xval*ones(size(Z));
        Y = yval*ones(size(Z));
        UU = [];
        VV = Z;
    end
else % sample_dim must be 2
    if numel(xval) == 1
        [Y, Z] = meshgrid(yy, zz);
        X = xval*ones(size(Y));
        UU = Y; VV = Z;
        u_offset = y_offset; v_offset = z_offset;
    end
    if numel(yval) == 1
        [X, Z] = meshgrid(xx, zz);
        Y = yval*ones(size(X));
        UU = X; VV = Z;
        u_offset = x_offset; v_offset = z_offset;
    end
    if numel(zval) == 1
        [X, Y] = meshgrid(xx, yy);
        Z = zval*ones(size(X));
        UU = X; VV = Y;
        u_offset = x_offset; v_offset = y_offset;
    end
end % if sample_dim
%% STEP 2: interpolate data f at all sample points X, Y, Z

% number of samples, for convenience

num_samp = numel(X);

% find elements containing mesh points and barycentric coordinates
% of mesh points within elements. 
%  [X(:), Y(:), Z(:)] is an (n x 3) matrix whose rows contain mesh points

[elt_samp, bary_samp] = tsearchn([VX',VY',VZ'], EToV, [X(:), Y(:), Z(:)]);

% elt_samp is a column vector of n element numbers and bary_samp is a
% (num_samp x 4) matrix of barycentric coordinates

% Convert barycentrics to r, s, t coordinate system: H&W p.409
% rst_samp is an (num_samp x 3) matrix whose rows contains local element 
% coordinates of sample points

rst_samp = bary_samp*[[-1, -1, -1]; [1, -1, -1]; [-1, 1, -1]; [-1, -1, 1]];

% Vandermonde matrix for evaluating field at each sample point
% (from field modal values on each sample element). V_samp is (num_samp x Np)
% i.e. each row holds the interpolation vector for one sample point.
% Note that Vandermonde3D behaves very strangely if it gets rows instead
% of columns for its r, s, t arguments!

V_samp = Vandermonde3D(N, rst_samp(:,1), rst_samp(:,2), rst_samp(:,3));

% Interpolation matrix for evaluating field at each sample point
% (from field nodal values on each sample element): global invV maps 
% from nodal to modal values. 

interp_samp = V_samp*invV;

% Perform interpolation for each sample value, saving results in interp.
% Maybe it would be more efficient to .* columns of interp_samp with rows
% of f(:,elt_samp) and sum the results? Probably not worth bothering about!

interp = zeros(num_samp, 1);
for j=1:num_samp
    interp(j) = interp_samp(j,:)*f(:,elt_samp(j));
end

% reshape flat interpolation results into same shape as U or V

if ~isempty(UU)
    interp = reshape(interp, size(UU));
elseif ~isempty(VV)
    interp = reshape(interp, size(VV));
end

%% STEP 3: plot data

if sample_dim == 1
    figure;
    if ~isempty(UU)
        plot(UU, interp)
    else
        plot(interp, VV)
    end
end

if sample_dim == 2
    figure;
    surf(UU, VV, interp);
    view(0,90)
    colorbar
    % OPTION 1
    % colormap(flipud(gray(2048)).^3)
    % OPTION 2
    cmap = gray;
    cmap = flipud(cmap(1:64,:));
    colormap(cmap);
    % following removes grid
    shading interp
end

% 24 July 2020
% colormap(flipud(gray(2048)).^3) does not work well for kernel plots
% flipping isn't great either for kernel plots

% # display data
% 
% if style == 'heat':
%     # matrices are indexed top-to-bottom and left-to-right wheras
%     # Cartesian coordinate are bottom-to-top and left-to-right hence
%     # flipup which reverses the order of the rows
%     plt.imshow(np.flipud(interp), 
%                extent=(U.min()-u_offset, U.max()+u_offset, 
%                        V.min()-v_offset, V.max()+v_offset), **kwargs)
% elif style == 'contour':
%     # horrible hack so levels can either be an integer (number of
%     # levels) or an array (level values)
%     if 'levels' in kwargs and isinstance(kwargs['levels'], int):
%         del kwargs['levels']
%         plt.contour(U, V, interp, N, **kwargs)
%     else:
%         plt.contour(U, V, interp, **kwargs)
% elif style == 'gradient':
%     # note np.gradient works out derivatives between rows then columns, 
%     # hence dV, dU rather than dU, dV
%     dV, dU = np.gradient(interp)
%     if reverse_gradient == True:
%         dU = -dU
%         dV = -dV
%     if unit_gradient:
%         nUV = np.sqrt(dU**2+dV**2)
%         dU, dV = dU/nUV, dV/nUV
%     plt.quiver(U, V, dU, dV, angles='xy')
% 
% return U, V, interp
% 
