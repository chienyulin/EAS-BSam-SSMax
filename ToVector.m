function v = ToVector(im)
% takes MxNx3 picture and returns (MN)x3 vector
sz = size(im);
if length(sz)==3
    v = reshape(im, [prod(sz(1:2)) sz(3)]);
else
    v = reshape(im, [prod(sz(1:2)) 1]);
end