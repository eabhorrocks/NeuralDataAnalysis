%% Firing rate comparison of different cell types

%% add these to your path
% https://github.com/eabhorrocks/NeuralDataAnalysis
% https://github.com/eabhorrocks/GenericFunctions

%% data path and files

% download link: https://figshare.com/articles/dataset/Data_for_Horrocks_et_al_2024/26031226?file=47033908
% put in a directory and update dataDir path below.
dataDir = 'D:\V1Data\Data\basic_111022'; % update to your correct directory

sessionTags = {
    'M22027', '20220517';
    'M22029', '20220607';
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'}

inputSuffix = 'basic.mat';

%% load and process each session to get firing rates of cells

tic % should take 2 mins

for isession = 1:5;
    fprintf('processing session: %i \n', isession)
    % load session
inputFileName = [sessionTags{isession,1},'_', sessionTags{isession,2},'_', inputSuffix];
load(fullfile(dataDir,inputFileName))


%%%%%% Process trials, partition into stat/run trials, stim/blanks %%%%%%
trialsSpeed2D = trials.Speed2D;
trialsSpeed2D(1) = []; % first trial can sometimes be unreliably rendered so we discard
tsd = trialsSpeed2D;

% split trials by state according to wheel data (V1 dynamics paper main criteria)
run_idx = find(cellfun(@(x) prop(x>0.5)>=0.75 & mean(x)>3, {tsd.WheelSpeed})); 
stat_idx = find(cellfun(@(x) prop(x<3)>=0.75 & mean(x)<0.5, {tsd.WheelSpeed}));

% assign trials as 'stat' or 'run', (and any others as 'nan')
[tsd.runFlag] = deal([nan]);
[tsd(stat_idx).runFlag] = deal([0]);
[tsd(run_idx).runFlag] = deal([1]);

tsd = tsd(~isnan([tsd.runFlag])); % remove trials that dont have an assigned behavioural state

% do some basic struct field renaming to work with getBinnedSpikeCounts
% function
for itrial =1 :numel(tsd)
    tsd(itrial).start_time = tsd(itrial).PDstart;
    tsd(itrial).absVel = abs(tsd(itrial).VelX1);
end

% get stimulus and blank trials to compare
blankTrials = tsd([tsd.numDots1]==0);
stimTrials = tsd([tsd.numDots1]==573 & [tsd.Contrast1]==1); % only full contrast trials

%%%%%% get trial-based spike counts for stimuli and blanks %%%%%%

for iunit = 1:numel(units)
    units(iunit).spiketimes = units(iunit).spike_times;
end

% get spike counts for stimTrials
options.intervalStart = 0; % start at stim onset
options.intervalEnd = 1; % get spikes in a window of 1s (stim duration)
options.binSpacing=1; % single 1s bin

[anUnits, cond] = getBinnedSpikeCounts(stimTrials, units, {'absVel', 'runFlag'}, options);

% assign relevant field to original 'units' struct array
for iunit = 1:numel(units)
    units(iunit).stimSpikes = anUnits(iunit).allSpikes;
    units(iunit).stimFR = mean(cellfun(@mean, units(iunit).stimSpikes),1); % FR for stat and run
end
% units(1).stimSpikes ia a cell array of spike counts 
% (6 speeds, 2 behavioural states)

% get spike counts for blankTrials
[anUnits, cond] = getBinnedSpikeCounts(blankTrials, units, {'absVel', 'runFlag'}, options);

% assign relevant field to original 'units' struct array
for iunit = 1:numel(units)
    units(iunit).blankSpikes = anUnits(iunit).allSpikes;
    units(iunit).blankFR = cellfun(@mean, units(iunit).blankSpikes); % FR for stat and run
end

% units(1).blankSpikes is a cell array of spike counts
% (1 stimulus condition, 2 behavioural states)

%%%%%% save processed units for later analysis %%%%%%
session(isession).units = units; 

end % finish looping through sessions
timeToProcessSessions = toc

%% combine all the units together and apply basic goodness criteria

% at this point, you may want to do this on a session-wise basis for
% appropriate stats. in which case you would just loop through
% session(isession).units

allUnits = cat(1,session.units);
nUnits = numel(allUnits)

% basic criteria from IBL white paper
goodUnits = allUnits([allUnits.isi_viol]<=0.1...
    & [allUnits.amplitude_cutoff]<=0.1 & [allUnits.amplitude]>=50);

% optionally filter out low firing rate cells
goodUnits = goodUnits([goodUnits.firing_rate]>=1);
nGoodUnits = numel(goodUnits)

%% split units into narrow- and wide-spiking and plot mean waveform

waveform_duration_criteria = 0.45; % from cell-explorer paper
narrow_idx = find([goodUnits.duration]<=0.45);
wide_idx = find([goodUnits.duration]>0.45);

nNarrow = numel(narrow_idx)
nWide = numel(wide_idx)

narrowWaveforms = cat(1,goodUnits(narrow_idx).mean_waveform);
wideWaveforms = cat(1,goodUnits(wide_idx).mean_waveform);

figure, hold on
shadedErrorBar(1:size(narrowWaveforms,2), mean(narrowWaveforms,1), sem(narrowWaveforms,1), 'lineProps', 'c')
shadedErrorBar(1:size(wideWaveforms,2), mean(wideWaveforms,1), sem(wideWaveforms,1), 'lineProps', 'm')
defaultAxesProperties(gca,true)


%% compare firing rates of narrow- and wide-spiking units 


narrow_spontaneousRate = cat(1,goodUnits(narrow_idx).blankFR); % stat, run
wide_spontaneousRate = cat(1,goodUnits(wide_idx).blankFR); % stat, run

narrow_stimRate = cat(1,goodUnits(narrow_idx).stimFR); % stat, run
wide_stimRate = cat(1,goodUnits(wide_idx).stimFR); % stat, run

%% plot results

% ignore state
figure, 
subplot(121), hold on
histogram(wide_spontaneousRate, 'BinWidth', 1, 'FaceColor', 'm','Normalization', 'Probability')
histogram(narrow_spontaneousRate, 'BinWidth', 1, 'FaceColor', 'c', 'Normalization', 'Probability')
xlim([0 50])
ylabel('Count'), xlabel('FR (Hz)')
title('Blank Trials')
subplot(122), hold on
histogram(wide_stimRate, 'BinWidth', 1, 'FaceColor', 'm','Normalization', 'Probability')
histogram(narrow_stimRate, 'BinWidth', 1, 'FaceColor', 'c','Normalization', 'Probability')
xlim([0 50])
xlabel('FR (Hz)')
title('Stimulus Trials')
sgtitle('Ignore state')


% stationary only
istate=1;
figure, 
subplot(121), hold on
histogram(wide_spontaneousRate(:,istate), 'BinWidth', 1, 'FaceColor', 'm','Normalization', 'Probability')
histogram(narrow_spontaneousRate(:,istate), 'BinWidth', 1, 'FaceColor', 'c', 'Normalization', 'Probability')
xlim([0 50])
ylabel('Count'), xlabel('FR (Hz)')
title('Blank Trials')
subplot(122), hold on
histogram(wide_stimRate(:,istate), 'BinWidth', 1, 'FaceColor', 'm','Normalization', 'Probability')
histogram(narrow_stimRate(:,istate), 'BinWidth', 1, 'FaceColor', 'c','Normalization', 'Probability')
xlim([0 50])
xlabel('FR (Hz)')
title('Stimulus Trials')
sgtitle('Stationary state')

% locomotion only
istate=2;
figure, 
subplot(121), hold on
histogram(wide_spontaneousRate(:,istate), 'BinWidth', 1, 'FaceColor', 'm','Normalization', 'Probability')
histogram(narrow_spontaneousRate(:,istate), 'BinWidth', 1, 'FaceColor', 'c', 'Normalization', 'Probability')
xlim([0 50])
ylabel('Count'), xlabel('FR (Hz)')
title('Blank Trials')
subplot(122), hold on
histogram(wide_stimRate(:,istate), 'BinWidth', 1, 'FaceColor', 'm','Normalization', 'Probability')
histogram(narrow_stimRate(:,istate), 'BinWidth', 1, 'FaceColor', 'c','Normalization', 'Probability')
xlim([0 50])
xlabel('FR (Hz)')
title('Stimulus Trials')
sgtitle('Locomotion state')
