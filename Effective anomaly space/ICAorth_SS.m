% Since IDICA-DR requires an initialization algorithm to generate a specific set if initial condition,
% the following My FastICA implements the FastICA using ATGP-generated data sample vectors 
% (program is called My_ATGP) as its initial condition
function [ICs,W,x_whitened]=ICAorth_SS(HIM,X,M)

bnd=size(HIM,3);
xx=size(HIM,1);
yy=size(HIM,2);

x=reshape(HIM,xx*yy,bnd);
x=x'; 

L=size(x,1);
K=size(x,2);

% %====Sphereing =====
u=mean(X,2); %dimension of u is 1*L
x_hat=X-u*ones(1,K);
C=(X*X')/size(X,2)-u*u';
% C=(x_hat*x_hat')/size(X,2);
%===========
[V,D]=eig(C);
D(D<0)=0;
A=pinv(sqrt(D))*V'; 
% A=inv(sqrt(D))*V';    % A is the whitening matrix....
 
x_whitened=A*(x_hat);
%  x_whitened=A*X;
%  x_whitened=X;
%====for cuprite data===
clear x;
clear x_hat;
%=========for initialization
[Loc,Sig]=My_ATGP(reshape(x_whitened',xx,yy,L),M);
W_initial=Sig; 
threshold = 0.0001;
B=zeros(L);
%===============find the first point=========

for round=1:M
    fprintf('IC %d', round); 
    %===Initial condition 
    w=W_initial(:,round);
    %===
    w=w-B*B'*w;
    w=w/norm(w);
    wOld=zeros(size(w));
    wOld2=zeros(size(w));
    i=1;
    while i<=3000
        w=w-B*B'*w;      
        w=w/norm(w);
%         fprintf('.');
        if norm(w-wOld)<threshold || norm(w+wOld)<threshold
            fprintf('Convergence after %d steps\n', i);
            B(:,round)=w;
            W(round,:)=w';
            break;
        end
        wOld2=wOld;
        wOld=w;
        w=(x_whitened*((x_whitened'*w).^3))/K-3*w;
        w=w/norm(w);
        i=i+1;
    end
    if (i>3000)
            fprintf('Warning! can not converge after 1000 steps \n, no more components'); 
            break;
    end
    round=round+1;    
end
ICs=W*x_whitened;
% W_hat=W*sqrt(D)*V';
% figure()
% for k=1:M 
%     s=reshape(abs(ICs(k,:)),xx,yy);
%     s=255*(s-min(min(s))*ones(size(s,1),size(s,2)))/(max(max(s))-min(min(s)));
%     temp=mean(reshape(s,xx*yy,1));
%     subplot(5,6,k); imshow(uint8(s));
% end
% A=W;
% x_re_ICA=V*sqrtm(D)*(A*ICs)+u*ones(1,xx*yy);
