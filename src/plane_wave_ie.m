function field_plane_wave = plane_wave_ie(t, x, y, z, param)

% this computes analytic plane wave propagating parallel to the vector
% [param.pwx, param.pwy, param.pwz]
% in the elastic case. This is inflated into the poroelastic data
% structure by adding zeros rows at positions 7 and 11:13

% e_vec_plane_wave_e is param.Nfields (13 or 9) x param.Nfields_e (9)
% c_plane_wave_ela has length param.Nfields_e (9)

e_vec_plane_wave_e = param.e_vec_plane_wave_e;
c_plane_wave_ela = param.c_plane_wave_ela;

% the expression for wavelength needs to be very carefully derived
wave_length = 2*pi*c_plane_wave_ela/param.omega_sou; 

% print eigenvalues to shell
if t==0
    fprintf(' PI phase speed: %4.2f\n', c_plane_wave_ela(9))
    fprintf(' S phase speed: %4.2f\n', c_plane_wave_ela(8))
    fprintf(' S phase speed: %4.2f\n', c_plane_wave_ela(7))
    fprintf(' PI wave length: %4.2f\n', wave_length(9))
    fprintf(' S wave length: %4.2f\n', wave_length(8))
    fprintf(' S wave length: %4.2f\n', wave_length(7))
end

%keyboard

% assemble analytic solutions
r1 = e_vec_plane_wave_e(:,9);%fast P
r2 = e_vec_plane_wave_e(:,8);%S-wave
r3 = e_vec_plane_wave_e(:,7);%S-wave

% component of x, y, z in direction of plane wave vector
p = (x*param.pwx+y*param.pwy+z*param.pwz) / ...
    sqrt(param.pwx^2+param.pwy^2+param.pwz^2);

% fast P-wave
tmp = exp(1i*param.omega_sou*(p/c_plane_wave_ela(9) - t));
field_plane_wave = zeros(size(tmp,1), size(tmp,2), param.Nfields);
field_plane_wave = addcomponent(tmp, r1, field_plane_wave, param.Nfields, param.weight(1));

% S-wave
tmp = exp(1i*param.omega_sou*(p/c_plane_wave_ela(8) - t));
field_plane_wave = addcomponent(tmp, r2, field_plane_wave, param.Nfields, param.weight(2));

% second S-wave
tmp = exp(1i*param.omega_sou*(p/c_plane_wave_ela(7) - t));
field_plane_wave = addcomponent(tmp, r3, field_plane_wave, param.Nfields, param.weight(3));

%keyboard

return

function field_plane_wave = addcomponent(tmp, r, field_plane_wave, Nfields, w)

for ii = 1:Nfields
    field_plane_wave(:,:,ii) = field_plane_wave(:,:,ii) + w*real(r(ii)*tmp);
end

return
