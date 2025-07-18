%% specify sessions

dataDir = 'D:\V1Data\Data\basic_111022';
sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

%% load session

isession = 1;
inputSuffix = 'basic.mat';
tic
inputFileName = [sessionTags{isession,1},'_', sessionTags{isession,2},'_', inputSuffix];
load(fullfile(dataDir,inputFileName))



%% get stat/run trials

trialsSpeed2D = trials.Speed2D;
trialsSpeed2D(1) = [];
tsd = trialsSpeed2D;

%split trials by state according to wheel data
run_idx = find(cellfun(@(x) prop(x>0.5)>=0.75 & mean(x)>3, {tsd.WheelSpeed}));
stat_idx = find(cellfun(@(x) prop(x<3)>=0.75 & mean(x)<0.5, {tsd.WheelSpeed}));

[tsd.runFlag] = deal([nan]);

[tsd(stat_idx).runFlag] = deal([0]);
[tsd(run_idx).runFlag] = deal([1]);

temp_tsd = tsd([tsd.numDots1]==573); % remove blnk trials.
temp_tsd = temp_tsd(~isnan([temp_tsd.runFlag]));
allStimConds = [vertcat(temp_tsd.VelX1), vertcat(temp_tsd.Contrast1), vertcat(temp_tsd.runFlag)];
[uniqueConds, ~, ic] = unique(allStimConds, 'rows');
tally = accumarray(ic, 1);
Result = [uniqueConds tally];

temp_tsd = temp_tsd([temp_tsd.Contrast1]==1); % only full contrast trials


%% get trial-based spike counts

for itrial =1 :numel(temp_tsd)
    temp_tsd(itrial).start_time = temp_tsd(itrial).PDstart;
    temp_tsd(itrial).absVel = abs(temp_tsd(itrial).VelX1);
end

for iunit = 1:numel(units)
    units(iunit).spiketimes = units(iunit).spike_times;
end


statTsd = temp_tsd([temp_tsd.runFlag]==0);
runTsd = temp_tsd([temp_tsd.runFlag]==1);
options.intervalStart = 0;
options.intervalEnd = 1.5;
options.binSpacing=1.5;

[anUnits, cond] = getBinnedSpikeCounts(temp_tsd, units, {'absVel', 'runFlag'}, options);

for iunit = 1:numel(units)
    units(iunit).allSpikes = anUnits(iunit).allSpikes;
end

for iunit = 1:numel(units)
    units(iunit).tuning = cellfun(@mean, units(iunit).allSpikes);
end

%% downsample to min # trials available
% set so equal # trials per condition

minTrial = min(min(cellfun(@(x) size(x,2), units(1).allSpikes)));

% units = units;

for iunit = 1:numel(units)
    units(iunit).allSpikesDownsample = cellfun(@(x) x(:,1:minTrial), units(iunit).allSpikes,'UniformOutput', false);
end

%% get cross-validated R2 (tuning strength) 
% options
r2opts.nPerms = 10;
r2opts.randFlag = 1;
r2opts.validMeans = 1;
r2opts.kval = 3;
r2opts.nShuffle = 100;

%get the tuning strength and significance
for iunit = 1:numel(units)

    %stationary trials
    sca = units(iunit).allSpikesDownsample(:,1)';
    [units(iunit).statR2, units(iunit).statR2_pval] = calc_kfold_R2(sca, r2opts.kval, r2opts.nPerms,...
        r2opts.randFlag, r2opts.validMeans, r2opts.nShuffle);

    %locomotion trials
    sca = units(iunit).allSpikesDownsample(:,2)';
    [units(iunit).runR2, units(iunit).runR2_pval] = calc_kfold_R2(sca, r2opts.kval, r2opts.nPerms,...
        r2opts.randFlag, r2opts.validMeans, r2opts.nShuffle);
end

%% get fano factor and dynamic range

for iunit = 1:numel(units)
    istate = 1;

    units(iunit).dynamicRange_stat = range(units(iunit).tuning(:,istate))/options.binSpacing;
    units(iunit).fanoFactor_stat = mean(cellfun(@(x) var(x)/mean(x), units(iunit).allSpikes(:,istate)),1);

    istate = 2;

    units(iunit).dynamicRange_run = range(units(iunit).tuning(:,istate))/options.binSpacing;
    units(iunit).fanoFactor_run = mean(cellfun(@(x) var(x)/mean(x), units(iunit).allSpikes(:,istate)),1);

end

%% fit descriptive functions

%clustering analysis reveals 4 main classes of tuning function

%gaussian fit and preferred speed with all trials

%stationary trials
istate = 1;
for iunit = 1:numel(units)

    %fit gaussians and find best fit
    [units(iunit).gaussParams_stat(1,:), units(iunit).gaussChar_stat, units(iunit).gaussR2_stat] = ...
        fitGaussianTemplates_tuning(units(iunit).allSpikes(:,istate),0.5,false);

    %get preferred stimulus
    if units(iunit).gaussChar_stat == 4 % if inverted, pref speed is min fr.
        [~, units(iunit).prefSpeed_stat] = min(units(iunit).tuning(:,istate));
    else % max fr
        [~, units(iunit).prefSpeed_stat] = max(units(iunit).tuning(:,istate));
    end

end

%locomotion trials
istate = 2;
for iunit = 1:numel(units)

    %fit gaussians and find best fit
    [units(iunit).gaussParams_run(1,:), units(iunit).gaussChar_run, units(iunit).gaussR2_run] = ...
        fitGaussianTemplates_tuning(units(iunit).allSpikes(:,istate),0.5,false);

    %get preferred stimulus
    if units(iunit).gaussChar_run == 4 % if inverted, pref speed is min fr.
        [~, units(iunit).prefSpeed_run] = min(units(iunit).tuning(:,istate));
    else % max fr
        [~, units(iunit).prefSpeed_run] = max(units(iunit).tuning(:,istate));
    end

end

%% Mutual information analysis

%MI options
optimiseBins = false; nMCSamples = 1000; sigFlag = false; correction = 'MLE'; nBinLim = 6; % nspeeds

for iunit = 1:numel(units)

    %stationary trials
    sca = units(iunit).allSpikesDownsample(:,1);

    [units(iunit).MI_stat, ~, ~, units(iunit).SSI_stat] =...
        calcMI(sca, correction, optimiseBins, sigFlag, nMCSamples, nBinLim);

    %locomotion trials
    sca = units(iunit).allSpikesDownsample(:,2);
    [units(iunit).MI_run, ~, ~, units(iunit).SSI_run] =...
        calcMI(sca, correction, optimiseBins, sigFlag, nMCSamples, nBinLim);

end

%% plots

r2_thresh = 0.1;
r2p_thresh = 0.05;

validIdx_stat = find(cat(1,units.statR2)>r2_thresh & cat(1,units.statR2_pval)<r2p_thresh);
validIdx_run = find(cat(1,units.runR2)>r2_thresh & cat(1,units.runR2_pval)<r2p_thresh);
validIdx_both = validIdx_stat(ismember(validIdx_stat, validIdx_run));

allTuning = cat(3,units.tuning)/options.binSpacing; % convert to spikes/s

% plot average tuning
figure
shadedErrorBar(1:6, mean(allTuning(:,1,validIdx_stat),3),...
    sem(allTuning(:,1,validIdx_stat),3),'lineProps','k');
hold on
shadedErrorBar(1:6, mean(allTuning(:,2,validIdx_run),3),...
    sem(allTuning(:,1,validIdx_run),3),'lineProps','r');
xlabel('Stimulus')
ylabel('Firing Rate (hz)')
title('Mean tuning curves')

figure
plot([units.statR2], [units.runR2],'k.')
hold on
plot([0.1 0.1], [-1 1],'r')
plot([-1, 1], [0.1 0.1],'r')
xlim([-0.6 1])
ylim([-0.6 1])
xlabel('Tuning Strength (stationary)')
ylabel('Tuning Strength (locomotion)')


% MI plots
figure
plot([units.MI_stat], [units.MI_run],'k.')
axis equal
hold on
plot([0 1.5], [0 1.5], 'r')
xlabel('Mutual information (stationary)')
ylabel('Mutual information (locomotion)')

figure
shadedErrorBar(1:6, mean(cat(1,units.SSI_stat)), sem(cat(1,units.SSI_stat)))
hold on
shadedErrorBar(1:6, mean(cat(1,units.SSI_run)), sem(cat(1,units.SSI_run)),'lineProps','r')
ylabel('SSI')
xlabel('Stimulus')

% preferred speeds joint histogram
allPrefs = [cat(1,units(validIdx_both).prefSpeed_stat), cat(1,units(validIdx_both).prefSpeed_run)];
vals = histcounts2(allPrefs(:,1), allPrefs(:,2));
vals = vals./sum(vals,2);
figure
imagesc(vals')
hold on
plot([0.5 6.5],[0.5 6.5],'r')
axis xy
colorbar, colormap(viridis)
xlabel('Preferred stimulus (stationary)')
ylabel('Preferred stimulus (locomotion)')

% plot gauss fit params for tuned bandpass cells
idx = find(cat(1,units.runR2)>r2_thresh & cat(1,units.runR2_pval)<r2p_thresh...
    & cat(1,units.statR2)>r2_thresh & cat(1,units.statR2_pval)<r2p_thresh...
    & cat(1,units.gaussChar_stat)==3 & cat(1,units.gaussChar_run)==3);

allStatParams = cat(1,units(idx).gaussParams_stat);
allRunParams = cat(1,units(idx).gaussParams_run);

minVals = min(cat(1,allStatParams,allRunParams));
maxVals = max(cat(1,allStatParams,allRunParams));

paramNames = {'baseline','amplitude','mu','sigma'};

figure
for iparam = 1:4
    subplot(2,2,iparam)
    plot(allStatParams(:,iparam), allRunParams(:,iparam),'k.')
    hold on
    plot([minVals(iparam), maxVals(iparam)],[minVals(iparam), maxVals(iparam)],'r')
    title(paramNames{iparam})
end


% plot dynamic range and fano factor
figure
plot([units(validIdx_both).dynamicRange_stat], [units(validIdx_both).dynamicRange_run],'k.')
xlabel('Dynamic range (stationary)')
ylabel('Dynamic range (locomotion)')

figure
plot([units(validIdx_both).fanoFactor_stat], [units(validIdx_both).fanoFactor_run],'k.')
xlabel('Fano Factor (stationary)')
ylabel('Fano Factor (locomotion)')

totalRunTime = toc