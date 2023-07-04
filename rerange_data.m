function HIM_new = rerange_data(HIM)
[m,n,l] = size(HIM);

HIM_2d = ToVector(HIM);
[HIM_2d_new,~,~] = scale_func(HIM_2d);

HIM_new = reshape(HIM_2d_new,m,n,l);



