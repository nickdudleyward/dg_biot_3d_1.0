function field_plane_wave = plane_wave_in(t, x, y, z, param)
% this computes analytic plane wave propagating parallel to the vector
% [param.pwx, param.pwy, param.pwz]

e_vec_plane_wave_in = param.e_vec_plane_wave_in;
c_plane_wave_in = param.c_plane_wave_in;

% the expression for wavelength needs to be very carefully derived
wave_length = 2*pi*c_plane_wave_in/param.omega_sou; 

% print eigenvalues to shell
if t==0
    fprintf(' PI phase speed: %4.2f\n', c_plane_wave_in(13))
    fprintf(' S phase speed: %4.2f\n', c_plane_wave_in(12))
    fprintf(' S phase speed: %4.2f\n', c_plane_wave_in(11))
    fprintf(' PII phase speed: %4.2f\n', c_plane_wave_in(10))
    fprintf(' PI wave length: %4.2f\n', wave_length(13))
    fprintf(' S wave length: %4.2f\n', wave_length(12))
    fprintf(' S wave length: %4.2f\n', wave_length(11))
    fprintf(' PII wave length: %4.2f\n', wave_length(10))
end

%keyboard

% assemble analytic solutions
r1 = e_vec_plane_wave_in(:,13);%fast P
r2 = e_vec_plane_wave_in(:,12);%S-wave
r3 = e_vec_plane_wave_in(:,11);%S-wave
r4 = e_vec_plane_wave_in(:,10);%slow P

% component of x, y, z in direction of plane wave vector
p = (x*param.pwx+y*param.pwy+z*param.pwz) / ...
    sqrt(param.pwx^2+param.pwy^2+param.pwz^2);

% fast P-wave
tmp = exp(1i*param.omega_sou*(p/c_plane_wave_in(13) - t));
field_plane_wave = zeros(size(tmp,1), size(tmp,2), param.Nfields);
field_plane_wave = addcomponent(tmp, r1, field_plane_wave, param.Nfields, param.weight(1));

% S-wave
tmp = exp(1i*param.omega_sou*(p/c_plane_wave_in(12) - t));
field_plane_wave = addcomponent(tmp, r2, field_plane_wave, param.Nfields, param.weight(2));

% second S-wave
tmp = exp(1i*param.omega_sou*(p/c_plane_wave_in(11) - t));
field_plane_wave = addcomponent(tmp, r3, field_plane_wave, param.Nfields, param.weight(3));

% slow P-wave
tmp = exp(1i*param.omega_sou*(p/c_plane_wave_in(10) - t));
field_plane_wave = addcomponent(tmp, r4, field_plane_wave, param.Nfields, param.weight(4));

%keyboard

return

function field_plane_wave = addcomponent(tmp, r, field_plane_wave, Nfields, w)

for ii = 1:Nfields
    field_plane_wave(:,:,ii) = field_plane_wave(:,:,ii) + w*real(r(ii)*tmp);
end

return

