function EAS_cube = EAS(HIM,para_set)

p = para_set.n_bs;

% determine j
[m,n,k] = size(HIM);
img0 = ToVector(HIM)';% L*N
[j,~,~,~,~,~,~,~,~]= MX_SVD(p,img0,m,n);

EAS_cube_temp = ICSC(HIM,img0,j);
EAS_cube = reshape(EAS_cube_temp',m,n,k);