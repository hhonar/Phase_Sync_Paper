% figure;
% matchix = matchstates(Corrsw20,Corrsin);
% count = 0;
% for ii = 1:2
%     for jj = 1:3
%         count = count + 1;
%         subplot(2,3,count);h = tight_subplot(2, 3, [.001 .001],[.01 .001],[.05 .05]);
% 
%         gsplot(CorrCV20(:,:,matchix(count)));
%         axis square; axis ij 
%         set(gca, 'XTick', [], 'YTick', [], 'CLim', [-1 1])%, 'XColor', [1 1 1], 'YColor', [1 1 1])
%         c = get(gca, 'Children');
%         set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');        
%        % text(1.5,-2,sprintf('State %d\n', ii, Sorder(ii)), 'Fontsize', 14);
%     end
% end


function [MATCHIDX,out] = matchstates(grnd,testmat,method)
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
    % [maxval,idx] = max(AA');
    % count = 0;
    % disp('level 1 - passed')
    % for i=1:numel(idx)
    %     if sum(idx == idx(i)) > 1
    %         count = count + 1
    %         ind{count} = find(idx == idx(i));
    %         val{count} = idx(i);
    %         
    %         [~,b] = max(maxval(ind{1}));
    %         idx(ind{count}(b)) = val{count};
    %         for j = 2:numel(ind)
    %             disp('level 2 - passed')
    %             ind{count}
    %             disp('level 3 - passed')
    %             ind{count}(j)
    %             disp('level 4 - passed')
    %             AA(ind{count}(j),:)
    %             disp('level 5 - passed')
    %           [aa,bb] = sort(AA(ind{count}(j),:),'descend');
    %           [C,IA] = setdiff(bb,idx);
    %           idx(ind{count}(j)) = C(1);
    %         end
    %     else
    %         idx = idx;
    %         
    %     end
    % end
    % 
    
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
    
     
        
%        sqrt(pdist2(reshape(CORRCOEF{1}(:,:,1),[],1)', reshape(CORRCOEF{2}(:,:,4),[],1)'))
%MATCHIDX = [idx;1:length(idx)]';

end
