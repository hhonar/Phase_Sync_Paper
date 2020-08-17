function [MATCHIDX,out] = matchstates(grnd,testmat,method)
%
%  Description: this dependecy is written to take:
%               grnd: a 3D matrix where the first two dimensions represent the covariance/correlation matrix, and the 3rd dim corresponds to state. (ground)
%               testmat: a 3D matrix of brain states stacked (:,:,1), ... (:,:,N) where N is number of the state (states to be matched with the grnd)
%               method: how to calculate the similarity between the cov/corr matrices of grnd and testmat
%
%
%
%
%

switch method
    case 1
        nStates = size(grnd,3);
        for j = 1:1:nStates
           for i = 1:1:nStates
               tmp = corrcoef(testmat(:,:,i),grnd(:,:,j));
               AA(j,i) = tmp(1,length(tmp));
           end  
        end

    
    case 2
        nStates = size(grnd,3);
        for j = 1:nStates
            for i = 1:nStates
                tmp = sqrt(pdist2(reshape(testmat(:,:,i),[],1)', reshape(grnd(:,:,j),[],1)'));
                AA(j,i) = -tmp(1,length(tmp));   % For correlations we wanted to find the max, but for the distance we want it to be min but since we want to use the same written function based on max to be used, simply it translates to find the max of -distance
            end
        end
        
        
    case 3
        nStates = size(grnd,3);
        for j = 1:nStates
            for i = 1:nStates        
                tmp = sqrt(pdist2(log2(reshape(testmat(:,:,i),[],1)'), log2(reshape(grnd(:,:,j),[],1)')));
                AA(j,i) = -tmp(1,length(tmp));   % For correlations we wanted to find the max, but for the distance we want it to be min but since we want to use the same written function based on max to be used, simply it translates to find the max of -distance
            end
        end

    
end
    count = 0;
    tmp = AA;
    for i = 1:nStates
    count = count + 1;
    maxall = max(tmp,[],'all');
    [i,j,v]=find(tmp == maxall,1,'last');
    idx(count,:) = [i j];
    tmp(i,:) = NaN;
    tmp(:,j) = NaN;
    tmp;
    end

    out = sortrows(idx,1);
    MATCHIDX = out(:,2);
    

end
