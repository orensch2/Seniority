%{
The purpose of this script is to analyze the concentration effect of local
seniority structure
%}

function measure=measure_local_seniority(x,k,S,y_grid)


%x=1; % magnitude of loss transmitted
%y=.2; % size of edge (bottle-neck parameter)
%k=3; % # of neighbours
%S=3; % # of seniority layers 
%y_grid=[x/k:.01:x]; %edge size
counter=1;

%{
Generate seniority structures:
each structure is a length k vector with integers between 0 and k that add
up to k
where each line is a unique structure
%}
seniority_structures=zeros(0,S);
for l=1:S
    aux=colex(k,l,1,k);
    aux=[aux, zeros(size(aux,1),k-size(aux,2))];
    seniority_structures(size(seniority_structures,1)+1:...
                         size(seniority_structures,1)+size(aux,1),:)=...
                         aux;
end

%{
Generate a vector of size k
Each element in the vector is the loss transmitted to a neighbour
%}
RESULT_stdev = NaN(size(y_grid,2),size(seniority_structures,1));
RESULT_entro = NaN(size(y_grid,2),size(seniority_structures,1));
for y=y_grid
    result=NaN(size(seniority_structures,1),k);
for structure=1:size(seniority_structures,1)
    
    result_per_node=zeros(k,1);
    for s=1:S
        loss=max(x-sum(seniority_structures(structure,s+1:end))*y,0); %total loss transmitted to layer s
        result_per_layer=min(loss/seniority_structures(structure,s),y);
        result_per_node(sum(seniority_structures(structure,1:s-1))+1:sum(seniority_structures(structure,1:s-1))+seniority_structures(structure,s))=...
            result_per_layer*ones(1,seniority_structures(structure,s));
    end
    result_per_node(isnan(result_per_node))=0;
    result(structure,:)=result_per_node/sum(result_per_node);
end

stdev=std(result,1,2);
entro=result.*(-log(result));
entro(isnan(entro))=0;
entro=1-sum(entro,2)/log(k); %1-entropy normalized by log(k)

RESULT_stdev(counter,:)=stdev;
RESULT_entro(counter,:)=entro;
counter=counter+1;
end

measure(:,:,1)=RESULT_stdev;
measure(:,:,2)=RESULT_entro;

end
