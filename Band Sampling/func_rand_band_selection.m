function HIM_bs = func_rand_band_selection(HIM,para_set)
% 
% randomly select para_set.n_bs bands.
%
% INPUTS:
%   - HIM:      hyperspectral imagecube (rows by columns by bands)
%   - para_set: parameters setting including maximum iteration number,
%               stopping rule criteria, and an initial random seed.
% OUTPUTS:
%   - HIM_bs:   an imagecube (rows by columns by n_bs)
%

n_bs = para_set.n_bs;
loop_rng_seed = para_set.loop_rng_seed;
rng(loop_rng_seed)

[~,~,l] = size(HIM);

temp_ind = 1:l;

ind_rand = temp_ind(randperm(length(temp_ind)));

HIM_bs = HIM(:,:,ind_rand(1:n_bs));
