%% Hierarchical sorting and clustering of neural responses
addpath(genpath('C:\Users\edward.horrocks\Documents\Code\GenericFunctions'))
%% load data

load('C:\Users\edward.horrocks\Documents\Code\NeuralDataAnalysis\datasets\exampleAllenPSTHs.mat');
% psth Array is a nResp x nTimeBin array of responses (each row is a PSTH).
% There's 500 reliable responses from each area (4000 total)
% areaVector is a vector of tags specifying the area a response was
% recorded in (VISp, VISl, VISal, VISam, VISpm, LGd, LP)

%% normalise responses to between 0 and 1

% to compare shapes of PSTHs, it's useful to normalise their firing rate. 
% Of course if you use a correlation distance measure or other measure that
% does not care about absolute firing rates then this may not be necessary.
allPSTHnorm = normalize(psthArray,2,'range');


%% do hierarchical sorting using dynamic time warping (DTW)

nPSTH = size(allPSTHnorm,1); % number of PSTH
dm = nan(nPSTH); % initialise dissimilarity matrix

maxStretch = 10; % this is a stretch factor for dynamic time warping. 
% Essentially the number of times an element can be repeated.


% generate dissimilarity matrix,
% nested loop version ran fastest when testing (~400 secs on 4 workers on
% old cpu)
tic 
for ipsth1 = 1:nPSTH
    psth_temp = allPSTHnorm(ipsth1,:);
    parfor ipsth2 = ipsth1:nPSTH
        % choose whichever distance measure suits your purposes!
        dm_temp(:,ipsth2) = dtw(psth_temp,allPSTHnorm(ipsth2,:), maxStretch); 
    end
    dm(ipsth1,:) = dm_temp;
    dm_temp = [];
end
toc

% dissimiarlity matrix is symmetric so we copy the upper triangular values
% to lower triangular elements of matrix
full_dm = triu(dm)+triu(dm,1)'; 

% generate a dendrogram using the full dissimilarity matrix
Z_sub = linkage(full_dm, 'average');

% get optimal leaf ordering of dendrogram (minimise sum of paired distances)
tic
leafOrder = optimalleaforder(Z_sub,full_dm);
toc % ~ 400s

% apply new ordering to the dendrogram to get the optimal sorting order of responses
[denHandle, leafnode_idx, outpermNodeOrder] = dendrogram(Z_sub, nPSTH, 'reorder', leafOrder, 'orientation', 'right');

% outpermNodeOrder is the order of the leaf nodes
% leafnode_idx is the label for individual responses to each leaf node

newPSTHorder = [];

for ileaf = 1:numel(outpermNodeOrder)
    idx = find(leafnode_idx==outpermNodeOrder(ileaf));
    newPSTHorder = vertcat(newPSTHorder, idx);
end

%% compare original order to sorted order
figure
ax1 = subplot(1,2,1);
imagesc(allPSTHnorm)
% defaultAxesProperties(gca, true)
ax = gca; ax.XTick = 1:20:200;
ax1.Colormap = crameri('bilbao');
title('original ordering')

ax2 = subplot(1,2,2);
imagesc(allPSTHnorm(newPSTHorder,:))
% defaultAxesProperties(gca, true)
ax = gca; ax.XTick = 1:20:200;
ax2.Colormap = crameri('bilbao');
title('optimal ordering')


%% get sorted dissimilarity matrix to compare to original

dm_sorted=nan;
% this is probably a dumb way of doing this. ~30s
tic
for ipsth1 = 1:nPSTH
    for ipsth2 = 1:nPSTH
        dm_sorted(ipsth1,ipsth2) = full_dm(newPSTHorder(ipsth1), newPSTHorder(ipsth2));
    end
end
toc

figure,
subplot(121)
imagesc(full_dm)
colorbar
title('original dissimilarity matrix')
subplot(122)
imagesc(dm_sorted), colormap(hot);
title('sorted dissimilarity matrix')
colorbar


%% some basic clustering

subclust_idx = cluster(Z_sub,'MaxClust',200); 

nSubClust = max(subclust_idx)

figure
h = histogram(subclust_idx,nSubClust);
idx = find(h.Values>50);
nHighCountClusters = numel(idx) % number of high count clusters
medCount = median(h.Values) % median cluster count

%% loop through high count clusters
for i = 1:numel(idx)
idx2 = find(subclust_idx==idx(i));
sequences = {};
for iseq =1:numel(idx2)
    sequences{iseq} = allPSTHnorm(idx2(iseq),:);
end
hold off
%average = DBA(sequences);
medidx = medoidIndex(sequences, maxStretch);
average = sequences{medidx};
plot(allPSTHnorm(subclust_idx==idx(i),:)')
hold on
plot(average, 'k', 'LineWidth', 2) 
hold on
% plot(mean(allPSTHnorm(subclust_idx==idx(i),:)), 'k--', 'LineWidth', 2)
ylim([0 1])
pause

end


