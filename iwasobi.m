function [W,Winit,ISR,signals]= iwasobi(x,AR_order,rmax,eps0)
%
% implements algorithm WASOBI for blind source separation of
% AR sources in a fast way, allowing separation up to 100 sources
% in the running time of the order of tens of seconds.
%
% INPUT:    x .... input data matrix d x N
%           d .... dimension of the data
%           N .... length of the data
%           ARmax .. maximum AR order of the separated sources
%           rmax ... a constant that may help to stabilize the algorithm.
%                    it has the meaning of maximum magnitude of poles of the
%                    AR sources. The choice rmax=1 means that no stabilization 
%                    is applied. The choice rmax=0.99 may lead to more stable
%                    results.
%           eps0 ... machine dependent constant to control condition number
%                    of weight matrices
%
% OUTPUT: W       ...... estimated de-mixing matrix
%         Winit ........ initial estimate of the matrix obtained by UWAJD
%         ISR .......... estimated ISR matrix which represents approximate accuracy 
%                        of the separation provided that there is no additive 
%                        noise in the model.
%         signals....... separated signals
%
% Code by Petr Tichavsky, using inputs from Eran Doron
% This version does not use ffdiag anymore.
%
if nargin<4
   eps0=5.0e-7;
end
if nargin<3
   rmax=0.99;
end  
num_of_iterations = 3;
[d N]=size(x);
Xmean=mean(x,2);
x=x-Xmean*ones(1,N);  %%%%%%%%%  removing the sample mean
T=length(x(1,:))-AR_order;
C0=corr_est(x,T,AR_order);
for k=2:AR_order+1
    ik=d*(k-1);
    C0(:,ik+1:ik+d)=0.5*(C0(:,ik+1:ik+d)+C0(:,ik+1:ik+d)');
end      %%%%%%%%% symmetrization
[Winit Ms] = uwajd(C0,20); %%% compute initial separation
                               %%% using uniform weights
%conver
%t1 = cputime-time_start;
W=Winit;
for in = 1:num_of_iterations
    [H ARC]=weights(Ms,rmax,eps0);
    [W Ms]=wajd(C0,H,W,5);
end
ISR=CRLB4(ARC)/N;
%t1 = [t1 cputime-time_start];
signals=W*x+(W*Xmean)*ones(1,N);

end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of IWASOBI

function G=THinv5(phi,K,M,eps)
%
%%%% Implements fast (complexity O(M*K^2))
%%%% computation of the following piece of code:
%
%C=[];
%for im=1:M 
%  A=toeplitz(phi(1:K,im),phi(1:K,im)')+hankel(phi(1:K,im),phi(K:2*K-1,im)')+eps(im)*eye(K);
%  C=[C inv(A)];
%end  
%
% DEFAULT PARAMETERS: M=2; phi=randn(2*K-1,M); eps=randn(1,2);
%   SIZE of phi SHOULD BE (2*K-1,M).
%   SIZE of eps SHOULD BE (1,M).

phi(2*K,1:M)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
almold=2*phi(1,:)+eps;
C0=1./almold;
x1=zeros(K,M); x2=x1; x3=x1; x4=x1;
x1(1,:)=C0; x2(1,:)=C0; 
x3(1,:)=-C0.*phi(2,:); 
x4(1,:)=-2*C0.*phi(2,:); 
x4old=[];
lalold=2*phi(2,:)./almold;
for k=1:K-1
    f2o=phi(k+1:-1:2,:)+phi(k+1:2*k,:);
    alm=sum(f2o.*x4(1:k,:),1)+phi(1,:)+eps+phi(2*k+1,:);
    a0=zeros(1,M); 
    if k<K-1
       a0=phi(k+2,:);
    end   
    gam1=sum(f2o.*x1(1:k,:),1);
    gam3=sum(f2o.*x3(1:k,:),1)+a0+phi(k,:);
    x4(k+1,:)=ones(1,M);
    b1m=sum(([phi(2:k+1,:); a0]+[zeros(1,M); phi(1:k,:)]).*x4(1:k+1,:));
    b2m=sum(([a0; phi(k+1:-1:2,:)]+phi(k+2:2*k+2,:)).*x4(1:k+1,:));
    latemp=b2m./alm;
    b2m=latemp-lalold; lalold=latemp;
    bom=alm./almold;
    ok=ones(k+1,1);
    x2(1:k+1,:)=x4(1:k+1,:).*(ok*(1./alm));
    x1(1:k+1,:)=[x1(1:k,:); zeros(1,M)]-(ok*gam1).*x2(1:k+1,:);
    x3(1:k+1,:)=[x3(1:k,:); zeros(1,M)]-(ok*gam3).*x2(1:k+1,:);
    x4temp=x4(1:k,:);
    x4(1:k+1,:)=[zeros(1,M); x4(1:k,:)]+[x4(2:k,:); ones(1,M); zeros(1,M)]...
       -(ok*bom).*[x4old; ones(1,M); zeros(1,M)]...
       -(ok*b2m).*x4(1:k+1,:)-(ok*b1m).*x1(1:k+1,:)-(ok*x4(1,:)).*x3(1:k+1,:);
    x4old=x4temp;
    almold=alm;
end
MK=M*K;
G=zeros(K,MK);
G(:,1:K:MK)=x1; clast=zeros(K,M);
f1=[phi(2:K,:); zeros(1,M)]+[zeros(1,M); phi(1:K-1,:)];
f2=[zeros(1,M); phi(K:-1:2,:)]+[phi(K+1:2*K-1,:); zeros(1,M)];
for k=2:K
    ck=G(:,k-1:K:MK);
    G(:,k:K:MK)=[ck(2:K,:); zeros(1,M)]+[zeros(1,M);  ck(1:K-1,:)]...
          -clast-(ok*sum(f1.*ck)).*x1-(ok*sum(f2.*ck)).*x2-(ok*ck(1,:)).*x3...
          -(ok*ck(K,:)).*x4;
    clast=ck;
end 

end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   of THinv5

function [AR,sigmy]=armodel(R,rmax)
%
% to compute AR coefficients of the sources given covariance functions 
% but if the zeros have magnitude > rmax, the zeros are pushed back.
%
[M,d]=size(R);
AR = zeros(M,d);
for id=1:d
    AR(:,id)=[1; -toeplitz(R(1:M-1,id),R(1:M-1,id)')\R(2:M,id)];
    v=roots(AR(:,id)); %%% mimicks the matlab function "polystab"
%    v1(1,id)=max(abs(v));
    vs=0.5*(sign(abs(v)-1)+1);
    v=(1-vs).*v+vs./conj(v);
    vmax=max(abs(v));
%    v2(1,id)=max(abs(v));
    if vmax>rmax
       v=v*rmax/vmax;
    end   
    AR(:,id)=real(poly(v)'); %%% reconstructs back the covariance function
end 
Rs=ar2r(AR);
sigmy=R(1,:)./Rs(1,:);
% [v1; v2]
end %%%%%%%%%%%%%%%%%%%%%%%  of armodel

function [ r ] = ar2r( a )
%%%%%
%%%%% Computes covariance function of AR processes from 
%%%%% the autoregressive coefficients using an inverse Schur algorithm 
%%%%% and an inverse Levinson algorithm (for one column it is equivalent to  
%%%%%      "rlevinson.m" in matlab)
% 
  if (size(a,1)==1)
      a=a'; % chci to jako sloupce
  end
  
  [p m] = size(a);    % pocet vektoru koef.AR modelu
  alfa = a;
  K=zeros(p,m);
  p = p-1;
  for n=p:-1:1
      K(n,:) = -a(n+1,:);
      for k=1:n-1
          alfa(k+1,:) = (a(k+1,:)+K(n,:).*a(n-k+1,:))./(1-K(n,:).^2);
      end
      a=alfa;
  end
%  
  r = zeros(p+1,m);
  r(1,:) = 1./prod(1-K.^2);
  f = r;
  b=f;
  for k=1:p 
      for n=k:-1:1
          K_n = K(n,:);
          f(n,:)=f(n+1,:)+K_n.*b(k-n+1,:);
          b(k-n+1,:)=-K_n.*f(n+1,:)+(1-K_n.^2).*b(k-n+1,:);
      end
      b(k+1,:)=f(1,:);
      r(k+1,:) = f(1,:);
  end       
end %%%%%%%%%%%%%%%%%%%%%%%%%%%  of ar2r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R_est=corr_est(x,T,q)
%
NumOfSources = size(x,1);
R_est = zeros(NumOfSources,(q+1)*NumOfSources);

for index=1:q+1
    R_est(:,NumOfSources*(index-1) + (1:NumOfSources)) = 1/T*(x(:,1:T)*x(:,index:T+index-1)');
end    

end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of corr_est
%
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

function [H ARC]=weights(Ms,rmax,eps0)
%
[d,Ld]=size(Ms);
L=floor(Ld/d);
d2=d*(d-1)/2;
R=zeros(L,d);
for index=1:L
    id=(index-1)*d;
    R(index,:)=diag(Ms(:,id+1:id+d)).';  %%% columns of R will contain 
                           %%% covariance function of the separated components
end
%
[ARC,sigmy]=armodel(R,rmax);      %%% compute AR models of estimated components
%
AR3=zeros(2*L-1,d2); 
ll = 1;
for i=2:d
  for k=1:i-1
      AR3(:,ll) = conv(ARC(:,i),ARC(:,k));
      ll = ll+1;
%    AR3=[AR3 conv(AR(:,i),AR(:,k))];
  end  
end
phi=ar2r(AR3);     %%%%%%%%%% functions phi to evaluate CVinv
H=THinv5(phi,L,d2,eps0*phi(1,:));  %%%% to compute inversions of CV 
                                       %%%% It has dimension zeros(M,M*d2).
im=1; 
for i=2:d
  for k=1:i-1
     fact=1/(sigmy(1,i)*sigmy(1,k));
     imm=(im-1)*L;
     H(:,imm+1:imm+L)=H(:,imm+1:imm+L)*fact;
     im=im+1;
  end
end  

end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of weights

function ISR = CRLB4(ARC)
%
% CRLB4(ARC) generates the CRLB for gain matrix elements (in term 
% of ISR) for blind separation of K Gaussian autoregressive sources 
% whose AR coefficients (of the length M, where M-1 is the AR order)
% are stored as columns in matrix ARC.

[M K]=size(ARC);

Rs=ar2r(ARC);

sum_Rs_s=zeros(K,K);

for s=0:M-1
    for t=0:M-1
        sum_Rs_s=sum_Rs_s+(ARC(s+1,:).*ARC(t+1,:))'*Rs(abs(s-t)+1,:);
    end
end

denom=sum_Rs_s'.*sum_Rs_s+eye(K)-1;
ISR=sum_Rs_s'./denom.*(ones(K,1)*Rs(1,:))./(Rs(1,:)'*ones(1,K));
ISR(eye(K)==1)=0;

end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of CRLB4

function [W_est Ms]=uwajd(M,maxnumiter,W_est0)
%
% my approximate joint diagonalization with uniform weights
%
% Input: M .... the matrices to be diagonalized, stored as [M1 M2 ... ML]
%        West0 ... initial estimate of the demixing matrix, if available
% 
% Output: W_est .... estimated demixing matrix
%                    such that W_est * M_k * W_est' are roughly diagonal
%         Ms .... diagonalized matrices composed of W_est*M_k*W_est'
%         crit ... stores values of the diagonalization criterion at each 
%                  iteration
% 
[d Md]=size(M);
L=floor(Md/d);
Md=L*d;
iter=0;
eps=1e-7;
improve=10;
if nargin<3
   [H E]=eig(M(:,1:d));
   W_est=diag(1./sqrt(diag(E)))*H';
else
   W_est=W_est0;
end  
if nargin<2
   maxnumiter=20;
end   
Ms=M;  
Rs=zeros(d,L);
for k=1:L
      ini=(k-1)*d;
      M(:,ini+1:ini+d)=0.5*(M(:,ini+1:ini+d)+M(:,ini+1:ini+d)');
      Ms(:,ini+1:ini+d)=W_est*M(:,ini+1:ini+d)*W_est';
      Rs(:,k)=diag(Ms(:,ini+1:ini+d));
end 
crit=sum(Ms(:).^2)-sum(Rs(:).^2);  
while improve>eps && iter<maxnumiter
  b11=[]; b12=[]; b22=[]; c1=[]; c2=[];
  for id=2:d 
    Yim=Ms(1:id-1,id:d:Md);
    b22=[b22; sum(Rs(id,:).^2)*ones(id-1,1)];
    b12=[b12; (Rs(id,:)*Rs(1:id-1,:)')'];
    b11=[b11; sum(Rs(1:id-1,:).^2,2)];
    c2=[c2; (Rs(id,:)*Yim')'];
    c1=[c1; sum(Rs(1:id-1,:).*Yim,2)];
  end
  det0=b11.*b22-b12.^2; 
  d1=(c1.*b22-b12.*c2)./det0; 
  d2=(b11.*c2-b12.*c1)./det0;
%    value=norm([d1; d2])
  m=0; 
  A0=eye(d);
  for id=2:d
      A0(id,1:id-1)=d1(m+1:m+id-1,1)';
      A0(1:id-1,id)=d2(m+1:m+id-1,1);
      m=m+id-1;
  end  
  Ainv=inv(A0);  
  W_est=Ainv*W_est;
  Raux=W_est*M(:,1:d)*W_est';
  aux=1./sqrt(diag(Raux));
  W_est=diag(aux)*W_est;  % normalize the result
  for k=1:L
     ini=(k-1)*d;
     Ms(:,ini+1:ini+d) = W_est*M(:,ini+1:ini+d)*W_est';
     Rs(:,k)=diag(Ms(:,ini+1:ini+d));
  end
  critic=sum(Ms(:).^2)-sum(Rs(:).^2);
%   improve=abs(critic-crit(end));
%   crit=[crit critic];
  improve=abs(critic-crit);
  crit = critic;
  iter=iter+1;
end  
end %%%%%%%%%%%%%%%%%%%  of uwajd

function [W_est Ms]=wajd(M,H,W_est0,maxnumit)
%
% my approximate joint diagonalization with non-uniform weights
%
% Input: M .... the matrices to be diagonalized, stored as [M1 M2 ... ML]
%        H .... diagonal blocks of the weight matrix stored similarly
%                     as M, but there is dd2 blocks, each of the size L x L
%        West0 ... initial estimate of the demixing matrix, if available
%        maxnumit ... maximum number of iterations
% 
% Output: W_est .... estimated demixing matrix
%                    such that W_est * M_k * W_est' are roughly diagonal
%         Ms .... diagonalized matrices composed of W_est*M_k*W_est'
%         crit ... stores values of the diagonalization criterion at each 
%                  iteration
%
%
[d Md]=size(M);
L=floor(Md/d);
dd2=d*(d-1)/2;
Md=L*d;
if nargin<4
   maxnumit=100;
end   
if nargin<3
   [H E]=eig(M(:,1:d));
   W_est=diag(1./sqrt(diag(E)))*H';
else
   W_est=W_est0;
end  
Ms=M;  
Rs=zeros(d,L);
for k=1:L
      ini=(k-1)*d;
      M(:,ini+1:ini+d)=0.5*(M(:,ini+1:ini+d)+M(:,ini+1:ini+d)');
      Ms(:,ini+1:ini+d)=W_est*M(:,ini+1:ini+d)*W_est';
      Rs(:,k)=diag(Ms(:,ini+1:ini+d));
end 
    
for iter=1:maxnumit
 b11=zeros(dd2,1); b12=b11; b22=b11; c1=b11; c2=c1;
 m=0; 
 for id=2:d        
    for id2=1:id-1
        m=m+1; im=(m-1)*L;
        Wm=H(:,im+1:im+L);
        Yim=Ms(id,id2:d:Md);
        Rs_id = Rs(id,:);
        Rs_id2 = Rs(id2,:);
        Wlam1=Wm*Rs_id';
        Wlam2=Wm*Rs_id2';
        b11(m)=Rs_id2*Wlam2;
        b12(m)=Rs_id*Wlam2;
        b22(m)=Rs_id*Wlam1;
        c1(m)=Wlam2'*Yim';
        c2(m)=Wlam1'*Yim';
     end
  end
  det0=b11.*b22-b12.^2; 
  d1=(c1.*b22-b12.*c2)./det0; 
  d2=(b11.*c2-b12.*c1)./det0;
  m=0; 
  A0=eye(d);
  for id=2:d
      A0(id,1:id-1)=d1(m+1:m+id-1,1)';
      A0(1:id-1,id)=d2(m+1:m+id-1,1);
      m=m+id-1;
  end  
  Ainv=inv(A0);  
  W_est=Ainv*W_est;
  Raux=W_est*M(:,1:d)*W_est';
  aux=1./sqrt(diag(Raux));
  W_est=diag(aux)*W_est;  % normalize the result
  for k=1:L
     ini=(k-1)*d;
     Ms(:,ini+1:ini+d) = W_est*M(:,ini+1:ini+d)*W_est';
     Rs(:,k)=diag(Ms(:,ini+1:ini+d));
  end
end 
end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of wajd

