function [X_ica]=ICSC(HIM,img0,j)

[no_lines,no_rows,no_bands]=size(HIM);
[ICs,Pw,x_whitened]=ICAorth_SS(HIM,img0,j);
Z=Pw;
X_ica=Z'*Z*x_whitened;
% sparse X_ica
T=X_ica;
[Temp,idx]=sort(abs(T(:)),'descend');
card=j*no_lines*no_rows;
Num=no_lines*no_rows*no_bands;
X_ica(idx(card+1:Num))=0;