clear; clc
%{
This script runs an endogenous seniority model in which nodes select into
securities with different seniority. A computer is needed in order to 
explore the model, in that payoff matrices of players depend on fixed 
points computed for the shocks assumed. 

The script assumes a shock distribution, computes a fixed point for each 
state of nature, and then gets expected moments for equity, contagion etc.

After reaching a full payoff matrix for each player the script computes 
Nash equilibrium.
%}

%parameters for fixed point iteration
precision=1e-5;
convergence=1; %arbitrary number > precision

%{
***************************************************************************
***************************BALANCE SHEET***********************************
***************************************************************************
%}

N=3; % # nodes
y=1; % edge weight
S=2; % # seniority layers
c_=.9; %outside assets
b_=.8; %outside liabilities
c=c_*ones(N*S,1);
b=b_*ones(N*S,1); %assuming outsiders are senior still


%{
***************************************************************************
***************************SHOCK DISTRIBUTION******************************
***************************************************************************
%}
shock_size=c_*1;
%X=[0,0,0; shock_size,0,0; 0,shock_size,0; 0,0,shock_size];
X=[zeros(1,N); shock_size*eye(N)];
no_shock_prob=.99;
f_X=[no_shock_prob,(1-no_shock_prob)/N*ones(1,N)];


%{
***************************************************************************
***************************NETWORK TOPOLOGY********************************
***************************************************************************
%}
Delta=ones(N)-eye(N); %complete network

%{
***************************************************************************
***************************NETWORK TOPOLOGY********************************
***************************************************************************
%}
possible_Actions=1:S;
aux=repmat({uint8(possible_Actions)},N*(N-1),1);
ALLCOMB=allcomb(aux{1:end}); %all types of seniority structures

%{
TABLE is the payoff matrix for each state (sen. struct, row) 
and for each player (column).

We run a loop for all possible states (parfor)
%}
TABLE=NaN(size(ALLCOMB,1),N);
for comb=1:size(ALLCOMB,1)
Action=ALLCOMB(comb,:); %each loop is for one seniority structure
expected_r=ones(S*N,1); %returns on securities
expected_r_=2*expected_r; %auxiliary variable for fixed pt.
while nansum(abs(expected_r-expected_r_))>precision
expected_r_=expected_r;


%{
***************************************************************************
******************SET UP THE LIABILITY MATRIX******************************
***************************************************************************
%}
%this liability matrix takes into consideration the expected recovery rates
%on securties and thus returns as well.
P=zeros(N,N,S);
%{
side note: the setting up of P(:,:,s) is transposed, this is all fine
the reason for that is the line P=reshape(P, N, N*S, 1)'; bellow
%}

%{
P(i,1+mod(i+j-1,N),Action(j+(i-1)*(N-1)))
liability of node i to node 1+mod(i+j-1,N) (the following j node, e.g. 
if j=1 for i=N-1 then the result is N, but for j=2 the result is 1.

The liability is at layer Action(j+(i-1)*(N-1))

y*expected_r(i+(Action(j+(i-1)*(N-1))-1)*N)

y is the face value of liability, which is then multiplied by the interest
i is paying on that seniority level.


%}

for i=1:N
    for j=1:N-1
        P(i,1+mod(i+j-1,N),Action(j+(i-1)*(N-1)))=...
            y*expected_r(i+(Action(j+(i-1)*(N-1))-1)*N);
    end
end

%{
***************************************************************************
***************************GET FIXED POINT*********************************
***************************************************************************
%}

P=reshape(P, N, N*S, 1)'; %reshaped to accord with notation in the paper
p_bar=sum(P,2); %total liabilities at each sen. layer

%get Pi -- relative liability matrix
Pi=repmat(bsxfun(@rdivide, P', p_bar'),S,1);
Pi(isnan(Pi))=0;

%get Pi transposed (can't simply transpose Pi because it has a stacked dim)
Pi_transpose=Pi(1:N,1:end);
Pi_transpose=reshape(Pi_transpose,N,N,S);
Pi_transpose=permute(Pi_transpose,[2 1 3]);
Pi_transpose=reshape(Pi_transpose, N, N*S, 1);
Pi_transpose=repmat(Pi_transpose,S,1);

%see paper for clarification of this matrix 
M=cell2mat(arrayfun(@(s) [repmat(eye(N),1,s-1),...
    zeros(N,N*(S+1-s))]',1:S,'un',0))';

p_per_shock=NaN(length(p_bar),size(X,1)); %payment vector per shock
r_per_shock=NaN(length(p_bar),size(X,1)); %recovery rate per shock
w_per_shock=NaN(N,size(X,1)); %node equity per shock

for shock=1:size(X,1) %loop for each state of the world, getting a fixed pt
c=c_*ones(N*S,1);
x=repmat(X(shock,:)',S,1);
c=c-x; %shock outside assets
    
p=p_bar;
convergence=1;
while convergence>precision

p_=max(min(p_bar, (c-b-M*p_bar+Pi*p)),0);
convergence=sum(abs(p_-p)); %l-1 norm of difference b/w iterations
p=p_;

end

%find equity for the fixed pt.
aux_=Pi_transpose*p;
aux_=aux_(1:N);
aux=cell2mat(arrayfun(@(i) sum(p(i:N:N*S)),1:N,'un',0))';
p_bar_=cell2mat(arrayfun(@(i) sum(p_bar(i:N:N*S)),1:N,'un',0))';
w=c(1:N)-b(1:N)+aux_-aux; %equity for the fixed pt

w=max(0,w); %if equity is negative this means outsiders are not paid in full
w_per_shock(:,shock)=w; %store value

r=p./p_bar; %recovery rate per security
p_per_shock(:,shock)=p;
r_per_shock(:,shock)=r; %recovery rates of each security, stored for each shock

end

expected_r=(1./(f_X*r_per_shock'))'; %1 over expected recovery = return

end
expected_w=(f_X*w_per_shock')';
%{
interest_paid_by=NaN(N,1);
interest_paid_to=NaN(N,1);

for i=1:N
%interest payment of node i
interest_paid_by(i)=nanmean(expected_r(i:N:N*S))-1;
%interest payment to node i
    aux=0;
    for j=1:N-1
        %Action(j+(i-1)*(N-1)) %gives a number in {1,...,S}
        %Action(j+(i-1)*(N-1)) is the action that i takes on j, ie the
        %security of the claim on the node that comes j steps after him
        %1+mod(i+j-1,N) the number of the node coming j after i
        
        %this is the return on the security of player 1+mod(i+j-1,N)
        %that player i chose
        aux=aux+expected_r(N*(Action(j+(i-1)*(N-1))-1)+1+mod(i+j-1,N));
    end
    interest_paid_to(i)=aux/(N-1)-1;
    

end
%[interest_paid_to, interest_paid_by, interest_paid_to-interest_paid_by]
%}
%store the payoffs for that seniority structure
TABLE(comb,:)=expected_w; %+y*(interest_paid_to-interest_paid_by);

end

%{
This part of the code generates an individual payoff matrix
it is redundant

aux=repmat({uint8(possible_Actions)},(N-1)^2,1);
player_possible_Actions=allcomb(aux{1:end});


comb__=NaN(N,size(ALLCOMB,1)/size(player_possible_Actions,1),size(player_possible_Actions,1));
parfor i=1:N
    aux=ALLCOMB;
    aux(:,1+(i-1)*(N-1):N-1+(i-1)*(N-1))=[];
    comb_=NaN(S^(N-1),size(player_possible_Actions,1));    
    for others_action=1:size(player_possible_Actions,1)
        comb=0;
        counter=1;
        while ~isnan(comb)
            if comb~=0
                aux(comb,:)=NaN(1,(N-1)^2);
            end
            [~,comb]=ismember(player_possible_Actions(others_action,:),...
                aux,'rows');
            if comb==0
                comb=[];
            else
                comb_(counter,others_action)=comb;
            end

            counter=counter+1;
            
        end
    end
    %size(comb_)
    comb__(i,:,:)=comb_;
end

TABLE=NaN(S^(N-1),size(player_possible_Actions,1),N);
for i=1:N
    for others_action=1:size(player_possible_Actions,1)
        TABLE(:,others_action,i)=TABLE_(comb__(i,:,others_action),i);    
    end
end
%}

%{
This part creates an array with N+1 dimensions, one dimension for each
player's strategy (of size S^(N-1)) and the last dimension for the payoffs
of players. Entry (i1,i2,...,iN,j) is the payoff of player j when player
k's strategy is ik. ik is a serial number of the strategy.
%}

%{
first create an array of NaN. To do so we use repmat of the matrix
NaN(S^(N-1),N)
along N-1 more dimensions, then permute it so that players' payoff
dimension is last (N+1)th dimension
%}
array=NaN(S^(N-1),N);
array=repmat(array,[1,1,S^(N-1)*ones(1,N-1)]);
array=permute(array,[1,N+1,3:N,2]);

%may need a flip to indx
%{
func provides a map from comb=1:size(ALLCOMB,1) to an N vector of entries
with the serial number of the strategy. Note that the order in which
ALLCOMB is arranged is important for this to work propperly.
%}
func = @(u,comb) mod(floor((comb-1)./(S^(N-1)).^(u-1)),(S^(N-1)))+1;

array_=NaN(size(ALLCOMB,1),N);
for i=1:N 
for comb=1:size(ALLCOMB,1)
%indx=mat2cell(func(1:N,comb),1,ones(N,1))
indx=num2cell(func(sort(1:N,'descend'),comb));
array_(comb,:)=func(sort(1:N,'descend'),comb);
array(indx{1:end},i)=TABLE(comb,i);
end
end


%{
This part finds the optimal response of each player to other actions
%}
%create a matrix aux of all possible strategies by others
aux=1:S^(N-1);
aux=repmat({uint8(aux)},N-1,1);
aux=allcomb(aux{1:end});
aux=mat2cell(aux,ones(1,size(aux,1)),ones(1,N-1));

%aux_ is a temporary vector of payoffs for a given player given others'
%actions
aux_=NaN(S^(N-1),1);  
%{
aux__ is a temporary vector of other players stra. as a serial number from
1 to S^(N-1). startegy number k refers to the kth row from temp
temp=1:S;
temp=repmat({uint8(aux)},N-1,1);
temp=allcomb(aux{1:end});
temp=mat2cell(aux,ones(1,size(aux,1)),ones(1,N-1));

where an entry col in that row refers to the type of security this player
has on the player that is col after him, that is if it is player i's
strategy then
mod(i+ciel(col/(N-1))-1,N)+1
%}
aux__=NaN(1,size(aux,2));
%{
indx_ is the serial number of the strategies of all players in the right
order. This is then mapped to comb_ by searching array_.
%}
indx_=NaN(1,N);
%{
ind is a cell array, with entry i collecting the set of states that, once 
reached, are desirable by player i.

ind_ is a temporary vector saving the desired states from the set that
player i is asked to choose from.
%}
ind=cell(N,1);
for i=1:N
    ind_=NaN(size(aux,1),1);
    for comb=1:size(aux,1)
        aux_(1:S^(N-1))=array(aux{comb,1:i-1},1:S^(N-1),aux{comb,i:end},i);
        [~,indx]=max(aux_);
        aux__(1:size(aux,2))=aux{comb,1:end};
        indx_(1:N)=[aux{comb,1:i-1},indx,aux{comb,i:end}];
        [~,comb_]=ismember(indx_,...
            array_,'rows');
        
        ind_(comb)=comb_;
        %get from comb (aux) +1 indx to comb of ALLCOMB
    end
    ind{i}=ind_;
end


%{
This part computes equilibrium by intersecting the set of acceptable states
for each player. Given that players -i made their choices, player i can
choose from a set of states of size S^(N-1).

That is, a full strategy of each player is a response (a function) that
takes the strategies of -i  and produces the desired strategy of i. Thus
for each column in the payoff matrix, player i would have chosen a row.
This would result in a particular state (seniority structure), which player
i does not wish to deviate from...

ind is N by 1 cell array, each cell containg a vector of desired states
once reached. To reach equilibrium we simply intersect those sets.

The resulting equilibrium is a vector (after transformation from cell) of
Nash eq. states, with the number corresponding to the seniority structure
in the row of matrix ALLCOMB
%}

equilibrium=intersect_multiple(ind);
equilibrium=cell2mat(equilibrium(1));
ALLCOMB(equilibrium,:)
aux=TABLE(equilibrium,:); %payoffs of each player at each equilibrium
%%
%{
for row_=1:size(TABLE,1)
    if sum(max(aux)>TABLE(row_,:))==4
        row_
    end
end

%}

