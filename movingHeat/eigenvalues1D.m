clear

%stencil:
%u_{-N}=0
%\dot u_{-j}=(u_{-j-1}-2u_{-j}+u_{-j+1})/h^2, j=N-1..1
%\dot u_{0}=2(u_{0}-u_{-1})/h^2-2(v_{0}-u_{0})/h
%\dot v_{0}=2(u_{1}-v_{0})/h^2+2(u_{0}-v_{0})/h
%\dot v_{j}=(v_{j-1}-2u_{j}+u_{j+1})/h^2, j=N-1..1
%v_{N}=0

alpha=0.01
L=1;
N=100;
h=L/N;
%u-block
%main laplacian
Au=diag(ones(N-1,1),-1)-2*diag(ones(N,1))+diag(ones(N-1,1),1); 
Au(1,1)=-1;
%coupling, u part
Au(N,:)=0; Au(N,N-1)=2; Au(N,N)=-2-alpha*h;
Buv=zeros(N,N);
Buv(N,1)=2*alpha*h;

%v-block
%main laplacian
Av=diag(ones(N-1,1),-1)-2*diag(ones(N,1))+diag(ones(N-1,1),1); 
Av(N,N)=-1;
%coupling, v part
Av(1,:)=0; Av(1,2)=2; Av(1,1)=-2-alpha*h;
Bvu=zeros(N,N);
Bvu(1,N)=2*alpha*h;

%combine blocks
A=[Au Buv; Bvu Av];
[v,l]=eig(A);
