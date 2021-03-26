% Basic routine to convert a file from a relative to an absolute path
% If the path begins with the local file separator(/ or \) or the home
% directory symbol (~), leave it alone. Otherwise, prepend the current
% working directory.
% SPE 2019/08/03
% Needed because file existence and, it seems, h5 file access will look
% along the path instead of in the current directory.

function p = absname(f)

if strcmp(f(1), filesep) || strcmp(f(1), '~')
    p = f;
else
    p = fullfile(pwd, f);
end
