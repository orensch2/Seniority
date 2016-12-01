clear; clc

%{
The purpose of this script is to find the clearing vector of a financial
network with a general strtucture of seniority
%}

%Primitives
S=2;
N=4;
c_=5; %outside assets
b_=4; %outside liabilities
c=c_*ones(N*S,1);
b=b_*ones(N*S,1); %assuming outsiders are senior still
precision=1e-5;


%shock
x_=c_; %shock size
i=1; %shocked node
x=x_*repmat([1:N]'==i,S,1);
c=c-x; %shock outside assets

%TODO: generalize my scipt for creation of a random graph
%Alternatively, simply generate S random graphs

%Simple adjacency matrix as an example, 
y=.001;
P=y*repmat(ones(N)-eye(N),S,1);
p_bar=sum(P,2);

Pi=repmat(bsxfun(@rdivide, P', p_bar'),S,1);
M=cell2mat(arrayfun(@(s) [ones(N,s-1), zeros(N,N*S+1-s)]',1:S,'un',0))';

convergence=1;
p=p_bar;
while convergence>precision

p_=max(min(p_bar, (c-b-M*p_bar+Pi*p)),0);
convergence=sum(abs(p_-p));
p=p_;


end
aux=Pi*p;
aux=cell2mat(arrayfun(@(i) sum(aux(i:N:N*S)),1:N,'un',0))';
w=c(1:N)-b(1:N)+aux








