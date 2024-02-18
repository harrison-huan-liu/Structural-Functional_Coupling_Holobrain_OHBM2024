function Graph = Thresholding_Code(Graph)
% for i = 1:size(Graph,2)
%     W = Graph(i).W;
    W = Graph;

    %% Set the thresholding parameters
    T = 0.95; % keeping 30% links

    W(W<0) = 0;
    [thr] = density_threshold(W,T);    
    W( W < thr ) = 0; 
    
%     Graph(i).W = W;
%     for j = 1:264
%         Graph(i).W(j,j) = 0;
%     end
    Graph=W;
    
    if ~isempty(find(sum(W,1)==0))
        disp(['Thresholding is too large .....'])
    end
% end
    
function [thr] = density_threshold(W,t)

N = size(W,1);
K = nnz(triu(W));

element_triu = triu(W);

kt = t * ((N^2-N)/2);
kt = ceil(kt);

if kt >=K
    thr = 0;
else
    [val,ind] = sort(element_triu(:),'descend');
    thr = val(kt);    
end
    
    

end

end