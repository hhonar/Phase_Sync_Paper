function [idx,M,C,BIC,AIC,DBI] = mykmeansup(CVUrow,nS,maxK)

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Doing the k-means clustering on the concatenated vectorized upper
% trianglular part of correlation matrices, obtained from sliding window
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% 1: running k-means for nS number of clusters, replicating 20 times, and
% using sample (selecting nS observations from CVUrow at random) to choose
% initial cluster centroid positions
[~,C,sumd] = kmeans(CVUrow,nS,'MaxIter',150,'Start','sample','Replicates',200); % replicates default 200

% 2: running k-means for nS number of clusters, with starting points as the
% centeroids matrix C found in the previous step [idx: centroid index]
[idx,C] = kmeans(CVUrow,nS,'MaxIter',1000,'Start',C); %1000 default

for i = 1:nS
    M(:,:,i) = vmconv(C(i,:),'vec2mat');
end


% For calculating the BIC for k-means
RSS = sum(sumd);

BIC = RSS + (log(size(idx,1)))*nS*size(C,2);
AIC = RSS + nS*size(C,2);

E = evalclusters(CSWVect,'kmeans','DaviesBouldin','klist',[1:maxK]);
DBI = E.OptimalK;