function [field_plane_wave_hf, mem_hf] = plane_wave_hf(t, x, y, z, param)
% this computes analytic plane wave propagating parallel to the vector
% [param.pwx, param.pwy, param.pwz] on the whole mesh (regardless of
% param.hf_elt)

e_vec_plane_wave_hf = param.e_vec_plane_wave_hf;
c_plane_wave_hf = param.c_plane_wave_hf;
ph_vel_hf = param.ph_vel_hf;
% used by nested function addcomponent
num_all_fields = param.Nfields + param.num_mem_hf;

% the expression for wavelength needs to be very carefully derived
wave_length_hf = 2*pi*c_plane_wave_hf/param.omega_sou; 

% print eigenvalues to shell
if t==0
    fprintf(' PI complex velocity: %4.2f %4.2fi \n', real(c_plane_wave_hf(16)), imag(c_plane_wave_hf(16)))
    fprintf(' S complex velocity: %4.2f %4.2fi \n', real(c_plane_wave_hf(15)), imag(c_plane_wave_hf(15)))
    fprintf(' S complex velocity: %4.2f %4.2fi \n', real(c_plane_wave_hf(14)), imag(c_plane_wave_hf(14)))
    fprintf(' PII complex velocity: %4.2f %4.2fi \n', real(c_plane_wave_hf(13)), imag(c_plane_wave_hf(13)))
    fprintf(' PI phase speed: %4.2f\n', ph_vel_hf (16))
    fprintf(' S phase speed: %4.2f\n', ph_vel_hf(15))
    fprintf(' S phase speed: %4.2f\n', ph_vel_hf(14))
    fprintf(' PII phase speed: %4.2f\n', ph_vel_hf(13))
    fprintf(' PI wave length: %4.2f\n', wave_length_hf(16))
    fprintf(' S wave length: %4.2f\n', wave_length_hf(15))
    fprintf(' S wave length: %4.2f\n', wave_length_hf(14))
    fprintf(' PII wave length: %4.2f\n', wave_length_hf(13))
end

% assemble analytic solutions
r1 = e_vec_plane_wave_hf(:,16);%fast P
r2 = e_vec_plane_wave_hf(:,15);%S-wave
r3 = e_vec_plane_wave_hf(:,14);%S-wave
r4 = e_vec_plane_wave_hf(:,13);%slow P

% component of x, y, z in direction of plane wave vector
p = (x*param.pwx+y*param.pwy+z*param.pwz) / ...
    sqrt(param.pwx^2+param.pwy^2+param.pwz^2);

% fast P-wave
tmp = exp(1i*param.omega_sou*(p/c_plane_wave_hf(16) - t));
field_plane_wave_hf = zeros(size(tmp,1), size(tmp,2), num_all_fields);
field_plane_wave_hf = addcomponent(tmp, r1, field_plane_wave_hf, param.weight(1));

% S-wave
tmp = exp(1i*param.omega_sou*(p/c_plane_wave_hf(15) - t));
field_plane_wave_hf = addcomponent(tmp, r2, field_plane_wave_hf, param.weight(2));

% second S-wave
tmp = exp(1i*param.omega_sou*(p/c_plane_wave_hf(14) - t));
field_plane_wave_hf = addcomponent(tmp, r3, field_plane_wave_hf, param.weight(3));

% slow P-wave
tmp = exp(1i*param.omega_sou*(p/c_plane_wave_hf(13) - t));
field_plane_wave_hf = addcomponent(tmp, r4, field_plane_wave_hf, param.weight(4));

% extract memory variables
mem_hf = field_plane_wave_hf(:,:,14:16);
field_plane_wave_hf = field_plane_wave_hf(:,:,1:13);

    function field_plane_wave_hf = addcomponent(tmp, r, field_plane_wave_hf, w)

    for ii = 1:num_all_fields
        field_plane_wave_hf(:,:,ii) = field_plane_wave_hf(:,:,ii) + w*real(r(ii)*tmp);
    end
    end % function field_plane_wave_hf

end



