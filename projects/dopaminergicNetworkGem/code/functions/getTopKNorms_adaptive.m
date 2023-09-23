function [A_rxns, adjustedA_norms, A_norms] = getTopKNorms_adaptive(A,k,flags)
%each time we select a row of the matrix, we want to bias ourselves to
%"independent" reactions. Such a problem is hard, but we use a heuristic
%method that subtracts off the best 

%Input:
%A: the covariance matrix of of samples
%k: the number of reactions to select
%flags: a 1/0 array telling which reactions are exchange reactions (1:
%not exch, 0: exch) So we will not select any reaction with a 1.

%Output:
%A_rxns: a list of indices corresponding to the top k reactions
%A_norms: the corresponding norms of A_rxns

if nargin < 3
    flags = zeros(size(A,2),1);
end

indices = [];
adjustedA_norms = [];
for i=1:k
    A_modified = remove_indices(A,indices);
    
    A_norms = sqrt(sum(A_modified.^2,2));
    
    [A_norms_sorted, A_index] = sort(A_norms,1,'descend');
%     A_index_flagged = A_index(find(~flags));
    % A_norms_sorted_flagged = A_norms_sorted(find(~flags));
    j=1;
    while flags(A_index(j))==1
       j=j+1; 
    end
    A_rxn = A_index(j);
    % A_norms_sorted = A_norms_sorted_flagged(1);
    
    indices = [indices; A_rxn];
    adjustedA_norms = [adjustedA_norms; A_norms_sorted(j)];
    fprintf('Reaction %d/%d selected\n', i,k);
end

A_rxns = indices;
A_norms = sqrt(sum(A.^2,2));
A_norms = A_norms(indices);


end

function [A] = remove_indices(A,indices)

    if isempty(indices)
        return;
    end
    
    %extract out the rows;
    B = A(indices,:);

    for i=1:size(A,2)
        x = (B'\A(i,:)')';
        A(i,:) = A(i,:)-x*B;
    end
    
    
end
