function mat = flirtmat_read(fname)
%flirtmat_read: get FLIRT affine matrix from file as written by flirt -omat
% Example:
%  mat = flirtmat_read(fname)
% See also: flirtmat_write, flirtmat2worldmat, worldmat2flirtmat

% Copyright 2009 Ged Ridgway <ged.ridgway gmail.com>

elements = textread(fname, '%f', 'delimiter', ' '); % read along rows
mat = reshape (elements,4,4); % stacked column wise, so needs:
mat = mat'; % transpose to give usual convention ([0 0 0 1] along 4th row)
