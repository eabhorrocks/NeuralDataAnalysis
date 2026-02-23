%% plot a tuning curve

load('M22027_20220517_basic.mat')

%% get binned spike counts for each stimulus condition of interest
options.intervalStart = 0;
options.intervalEnd = 1.5; % set to trial duration
options.binSpacing=1.5; % set to trial duration
options.spikeTimesField = 'spike_times';

trialsStruct = trials.Speed2D;
paramsOfInterest = {'VelX1'};

[anUnits, cond] = getBinnedSpikeCounts(trials.Speed2D, units, {'VelX1'}, options);

for iunit = 1:numel(units)
    units(iunit).allSpikes = anUnits(iunit).allSpikes;
end

for iunit = 1:numel(units)
    units(iunit).tuning = cellfun(@mean, units(iunit).allSpikes);
end

%%
figure
errorbar(1:7, cellfun(@(x) mean(x,2), units(300).allSpikes), cellfun(@(x) sem(x,2), units(300).allSpikes))
ylabel('Spike count')
xlabel('Stimulus Value')

