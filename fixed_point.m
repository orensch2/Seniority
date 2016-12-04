clear; clc

%{
The purpose of this script is to find the clearing vector of a financial
network with a general strtucture of seniority
%}

%Primitives
N=4;
S=3; %has to be <=N-1
c_=1; %outside assets
b_=0; %outside liabilities
c=c_*ones(N*S,1);
b=b_*ones(N*S,1); %assuming outsiders are senior still
precision=1e-5;

%{
Not accounting for repetitions of structures...
seniority_structures=colex(N-1,S,0,N-1);
seniority_structures=seniority_structures(1:end-1,:); %get rid of last line, as (k,0) and (0,k) is same
%}
%{
Generate seniority structures:
each structure is a length k=N-1 vector with integers between 0 and k that add
up to k
where each line is a unique structure
%}
seniority_structures=zeros(0,N-1);
for l=1:S
    aux=colex(N-1,l,1,N-1);
    aux=[aux, zeros(size(aux,1),N-1-size(aux,2))];
    seniority_structures(size(seniority_structures,1)+1:...
                         size(seniority_structures,1)+size(aux,1),:)=...
                         aux;
end

seniority_structure=2; %choose 1 seniority structure (line) from matrix

%shock
x_=c_*3; %shock size
i=1; %shocked node
x=x_*repmat([1:N]'==i,S,1);
c=c-x; %shock outside assets

%TODO: generalize my scipt for creation of a random graph
%Alternatively, simply generate S random graphs

%Generate S liability matrices
%Each neighbour holds exactly one type of security (edge) on another
%All edges are of the same size y
P=zeros(N,N,S);
for i=1:N
    for s=1:S
        first=i+1+sum(seniority_structures(seniority_structure,1:s-1));
        first=1+mod(first-1,N);
        last=first-1+seniority_structures(seniority_structure,s);
        last=1+mod(last-1,N);
        if first~=i && last~=i
        if first<=last && ~(first==1 && last==N)
            P(i,first:last,s)=1;
        elseif first>last && (first~=last+1)
            P(i,first:N,s)=1;
            P(i,1:last,s)=1;
        end
        end
    end
end

y=1.5;
P=P*y;
P=reshape(P, N, N*S, 1)'; %reshaped to accord with notation in the paper
p_bar=sum(P,2); %total liabilities at each sen. layer

Pi=repmat(bsxfun(@rdivide, P', p_bar'),S,1);
Pi(isnan(Pi))=0;
Pi_transpose=Pi(1:N,1:end);

Pi_transpose=reshape(Pi_transpose,N,N,S);
Pi_transpose=permute(Pi_transpose,[2 1 3]);
Pi_transpose=reshape(Pi_transpose, N, N*S, 1);
Pi_transpose=repmat(Pi_transpose,S,1);

M=cell2mat(arrayfun(@(s) [repmat(eye(N),1,s-1), zeros(N,N*(S+1-s))]',1:S,'un',0))';

convergence=1; %arbitrary number > precision
p=p_bar;
while convergence>precision

p_=max(min(p_bar, (c-b-M*p_bar+Pi_transpose*p)),0);
convergence=sum(abs(p_-p)); %l-1 norm of difference b/w iterations
p=p_;


end
aux_=Pi_transpose*p;
aux_=aux_(1:N);
aux=cell2mat(arrayfun(@(i) sum(p(i:N:N*S)),1:N,'un',0))';
p_bar_=cell2mat(arrayfun(@(i) sum(p_bar(i:N:N*S)),1:N,'un',0))';
w=c(1:N)-b(1:N)+aux_-aux;

w_before_shock=(c_-b_)*ones(N,1);
[w_before_shock, w]
loss_dist=w_before_shock-w;




%[p_bar, p]
%[aux, aux_]

