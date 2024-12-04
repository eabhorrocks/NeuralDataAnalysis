%% Example dimensionality reduction using DataHigh toolbox

%% load dataset

load('Z:\ibn-vision\DATA\SUBJECTS\M24016\analysis\20240626\session_clusters_RUN1.mat')

% addpath to https://github.com/eabhorrocks/GenericFunctions

%% pre-process and filter data

session_clusters.spatial_response; % these ones are pre-smoothed firing rate, use unsmoothed for gpfa
% each row of cell array is a cluster
% each row of cell is a lap for that cluster
% each column of cell is a spatially binned response

cluster_regions = [session_clusters.region];
areas2use = {'V1_L','V1_R','HPC_R','HPC_L'};
minFR = 1;

lapCounts = cellfun(@(x) size(x,1), session_clusters.spatial_response(1,:));
trackIDs = repelem([1,2], 1, lapCounts);


% combine lap1 and lap2 to make zscoring easier later
for icluster = 1:size(session_clusters.spatial_response)
    clusters{icluster} = cat(1,session_clusters.spatial_response{icluster,1},session_clusters.spatial_response{icluster,2});
end
meanFRs = cellfun(@(x) mean(x,'all'), clusters); % mean spike count for each spatial bin

% z-score respnses to stop indiviudal units fro dominating latent factors
clusters_z = cellfun(@(x) zscore(x,[],'all'), clusters, 'UniformOutput', false);

% use clusters with mean spike count > minFR and from chosen region
clusters2use_idx = find(meanFRs(:)>=minFR & ismember([session_clusters.region],areas2use));
clusters2use = clusters_z(clusters2use_idx);


D=struct();
for ilap = 1:sum(lapCounts) % loop through laps for each track
    for icluster = 1:numel(clusters2use) % loop through clusters (for now, loop)
        temp = clusters2use{icluster};
        D(ilap).data(icluster,:) = temp(ilap,:);
        D(ilap).condition = trackIDs(ilap);
    end
end



%% Dim reduction
handles = []; % no longer needed
alg = 3; % 1PCA, 2PPCA, 3FA, 4 LDA, 5 GPFA

% cross-validate to find dimensions that maximise likelihood
candidateDims = 5:50; % note, GPFA cross-validaiton is slow, so be smart about the dim search

[projs, mse, like] = cvreducedims_edd(D, alg, candidateDims, handles);
[~, idx] = max(like); % find q that maximises likelihood of data
q = candidateDims(idx);

% fit FA with cross-validated dims
[newD, C, lat, explained, params] = reducedims_EH(D,alg, q, handles); % do FA using q dims

% info about dim reduction model
s.q = q;
s.qOpt = compute_dshared(params);
s.LoadingSim = compute_load_sim(params);
s.SV = compute_perc_shared(params);
s.params = params;
s.propSharedVariance = [lat(1); diff(lat)];
s.loadings = C;
s.nUnits = size(C,1);
s.D=D;
s.explained = explained;

% format trajectories into useful format
for itrial = 1:numel(newD)
    newD(itrial).y = D(itrial).data;
end

for itrack = 1:numel(lapCounts)
    idx = find([newD.condition]==itrack);
    cond(itrack).catData = cat(3,newD(idx).data);
    cond(itrack).catData = permute(cond(itrack).catData, [2, 1, 3]);
    cond(itrack).meanTrajectory = mean(cond(itrack).catData,3);
    cond(itrack).semTrajectory = sem(cond(itrack).catData,3);
end

s.cond= cond;

%% plot some basic metrics

% plot shared variance
figure
subplot(131)
plot(s.propSharedVariance,'k', 'Marker','.'), grid on
xlim([0.5 q+0.5]);
ylabel('Proportion of shared variance explained')
xlabel('Factor rank')
title(['Proportion of variance shared: ', num2str(s.SV)])

subplot(132)
imagesc(s.loadings);
ylabel('neuron')
xlabel('latent factor')
colormap(redblue)
maxVal = max(abs(s.loadings(:)));
caxis([-maxVal, maxVal]), colorbar
title(['qOpt: ', num2str(s.qOpt)])

subplot(133)
plot(s.LoadingSim, 'k', 'Marker','.'), grid on
xlim([0.5 q+0.5]);
ylabel('Loading Similarity')
xlabel('Factor rank')

%% plot mean trajectories for top factors

% assumes 2 tracks here...if you have more, index into them as
% cond(itrack).meanTrajectory. Same goes for other trajectory plots

figure
for idim = 1:8
    subplot(2,4,idim), hold on
    shadedErrorBar(1:size(s.cond(1).meanTrajectory,1), s.cond(1).meanTrajectory(:,idim),...
        cond(1).semTrajectory(:,idim),'lineProps','r')
    shadedErrorBar(1:size(s.cond(2).meanTrajectory,1), s.cond(2).meanTrajectory(:,idim),...
        cond(2).semTrajectory(:,idim),'lineProps','b')
    title(['Dim: ', num2str(idim), ', SV: ', num2str(s.propSharedVariance(idim)*100,2),'%'])
end


%% plot some single trials for a specific dim

idim = 1;

figure, hold on
plot(squeeze(s.cond(1).catData(:,idim,:)),'r');
plot(squeeze(s.cond(2).catData(:,idim,:)),'b');
title(["Dim: ", num2str(idim)])


%% 3d plot 

dims2plot=1:3;

figure, hold on
plot3(s.cond(1).meanTrajectory(:,dims2plot(1)),s.cond(1).meanTrajectory(:,dims2plot(2)),s.cond(1).meanTrajectory(:,dims2plot(3)),'r')
plot3(s.cond(2).meanTrajectory(:,dims2plot(1)),s.cond(2).meanTrajectory(:,dims2plot(2)),s.cond(2).meanTrajectory(:,dims2plot(3)),'b')
view(-30,15), grid on
xlabel(dims2plot(1)), ylabel(dims2plot(2)), zlabel(dims2plot(3))
