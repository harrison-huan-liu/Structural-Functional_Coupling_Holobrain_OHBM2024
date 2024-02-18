function [HW,thr,knum] = density_thresholding(W,t,knum)

%%% W adjancent matrix
%%% t density
%%% knum keep the link number for each node

W(W<0) = 0;
N = size(W,1);
K = nnz(triu(W));

kt =  t * ((N^2-N)/2);
kt = ceil(kt);

if kt > K
    
    thr = 0;
    HW = W;
    knum = knum;
    
else
    
    if knum > kt/(1.5*N)   
        
        knum = ceil(kt/(1.5*N));
        
    end
    
    MSK = zeros(N,N);
    for i = 1 : N
        
        tmp = W(:,i);
        knz = nnz(tmp);
        if knz <= knum

            ind = find(tmp > 0);
            MSK(ind(:),i) = 1;

        else

            [~,index] = sort(tmp,'descend');
            ind = index(1:knum);
            MSK(ind(:),i) = 1;

        end

    end

    MSK = (MSK + MSK') / 2;
    MSK(MSK>0) = 1;
    
    Km = nnz(triu(MSK));
    
    Ktt = kt - Km;
    
    Wk = W .* MSK;
    Ws = W .* (1 - MSK);
    element_triu = triu(Ws);
    [val,~] = sort(element_triu(:),'descend');
    
    thr = val(Ktt);
    Ws(Ws<thr) = 0;
    HW = Ws + Wk;
    HW = (HW + HW')/2;
    
end

end
% 
% function [thr] = density_threshold(W,t)
% 
% N = size(W,1);
% K = nnz(triu(W));
% 
% element_triu = triu(W);
% 
% kt =  t * ((N^2-N)/2);
% kt = ceil(kt);
% 
% if kt >=K
%     thr = 0;
% else
%     [val,ind] = sort(element_triu(:),'descend');
%     thr = val(kt);    
% end
%     
%     
% 
% end