function g = ricker_wavelet(t, t0_sou, frequency)

% second derivative of Gaussian
c = (pi*frequency)^2;
g = (1 - 2*c*(t-t0_sou).^2 ).*exp(-c*(t-t0_sou).^2);
