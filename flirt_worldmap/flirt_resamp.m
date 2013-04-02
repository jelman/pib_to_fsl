function flirt_resamp(affmat, src, trg, outname)
%flirt_resamp: Resample a NIfTI source image using a flirt -omat matrix
% Example:
%  flirt_resamp('src_flirtmat.txt', 'src.nii', 'trg.nii', 'res.nii')
% Note that this was written for testing purposes -- it is very slow, and
% better results could be achieved using e.g. spm_bsplinc/spm_bsplins.
% See also: flirtmat2worldmat

% Copyright 2009 Ged Ridgway <ged.ridgway gmail.com>

src = nifti(src);
trg = nifti(trg);

[worldmat spmvoxmat] = flirtmat2worldmat(affmat, src, trg);

res = trg;
res.dat.fname = outname;
create(res);
res.dat(:, :, :) = resamp_vox(spmvoxmat, src.dat(:, :, :), trg.dat.dim);

%% 
function res = resamp_vox(mat, src, trgdim)
% mat should be vox-vox and zero-based
x = 1:trgdim(1);
y = 1:trgdim(2);
z = 1:trgdim(3);
[x y z] = ndgrid(x, y, z);
xyz1 = [x(:) y(:) z(:) ones(size(x(:)))]';
XYZ1 = mat * xyz1; % transform to source coords
X = reshape(XYZ1(1, :), trgdim);
Y = reshape(XYZ1(2, :), trgdim);
Z = reshape(XYZ1(3, :), trgdim);
res = interpn(src, X, Y, Z);
