function flirtmat_write(fname, mat)
%flirtmat_write: save a 4-by-4 matrix to a file as handled by flirt -init
% Example:
%  flirtmat_write(fname, affinemat)
% See also: flirtmat_read, flirtmat2worldmat, worldmat2flirtmat

% Copyright 2009 Ged Ridgway <ged.ridgway gmail.com>

fid = fopen(fname, 'w');
% Note that transpose is needed to go from MATLAB's column-major matrix to
% writing the file in a row-major way (see also flirtmat_read)
str = sprintf('%f  %f  %f  %f\n', mat');
fprintf(fid, str(1:end-1)); % drop final newline
fclose(fid);
