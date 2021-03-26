function s = ftime(t)

% Formatted an elapsed time in seconds into days, hours, minutes
% and seconds. Returns a string.

% names of time units and number of each in the next unit up
% 1e15 is a virtual infinity, intended to ensure that the loop
% always stops after the last unit

n = {'s', 'm', 'h', 'd'};
u = [60 60 24 1e15];

% build the formatted string from the right, the least significant unit,
% processing quotient and remainder until we hit a zero quotient
% (the 1e15 at the end of u is intended to force a zero quotient when we
% reach the end of the specified units)

s = '';
for j=1:numel(u)
    q = floor(t/u(j));
    r = t - q*u(j);
    s = sprintf('%g%s %s', r ,n{j}, s);
    t = q;
    if t == 0
        break
    end
end
s = s(1:end-1);


