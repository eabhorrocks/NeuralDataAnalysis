function [units, cond, varInfo] = getBinnedSpikeCounts(trials, units, varsOfInterest, options)
% getBinnedSpikeCounts: Bins spike times for multiple units across different trial conditions.
%
% [units, cond, varInfo] = getBinnedSpikeCounts(trials, units, varsOfInterest, options)
%
% This function takes raw spike time data and trial information and returns
% binned spike counts organized by experimental conditions. It modifies the
% input 'units' struct by adding a field with the binned spike counts and
% also returns a separate 'cond' struct with detailed trial information.
%
% Inputs:
%   trials - (struct array) Each element must contain fields corresponding to
%            the variables of interest and a trial start time.
%   units  - (struct array) Each element must contain a field with spike times.
%            The name of this field can be specified in options.spikeTimesField.
%   varsOfInterest - (cell array of strings) Names of the fields in 'trials'
%                    that define the experimental conditions (e.g., {'Dir', 'Speed'}).
%   options - (struct) Optional parameters:
%     .intervalStart    - (scalar) Start of binning window relative to trial start (default: 0).
%     .binSpacing       - (scalar) Width of each bin in seconds (default: 0.1).
%     .intervalEnd      - (scalar) End of binning window relative to trial start (default: 1).
%     .uniqueVals       - (cell array) Pre-specified unique values for each varOfInterest.
%     .trialStartField  - (string) Name of the field in 'trials' for trial start time (default: 'start_time').
%     .spikeTimesField  - (string) Name of the field in 'units' for spike times (default: 'spiketimes').
%     .allSpikesField   - (string) Name of the output field to be added to 'units' (default: 'allSpikes').
%
% Outputs:
%   units   - (struct array) The original 'units' struct with an added field
%             (name specified by .allSpikesField) containing the binned spike counts
%             in a multi-dimensional cell array.
%   cond    - (struct array) Contains detailed information for each condition,
%             including the spike counts for all units for that condition.
%   varInfo - (struct array) Information about the variables of interest, including
%             their names, unique values, and corresponding dimension in the output.
%
% Example:
%   varsOfInterest = {'Direction', 'Speed'};
%   options.allSpikesField = 'binned_responses';
%   [units, cond, varInfo] = getBinnedSpikeCounts(myTrials, myUnits, varsOfInterest, options);

% --- 1. Set up Options and Parameters ---
if ~exist('options','var'),             options=struct;                 end
if ~isfield(options,'intervalStart'),   options.intervalStart = 0;      end
if ~isfield(options,'binSpacing'),      options.binSpacing = 0.1;       end
if ~isfield(options,'intervalEnd'),     options.intervalEnd = 1;        end
if ~isfield(options,'uniqueVals'),      options.uniqueVals = [];        end
if ~isfield(options,'trialStartField'), options.trialStartField = 'start_time'; end
if ~isfield(options,'spikeTimesField'), options.spikeTimesField = 'spiketimes'; end
if ~isfield(options,'allSpikesField'),  options.allSpikesField = 'allSpikes';   end

nVars = numel(varsOfInterest);
nUnits = numel(units);

% --- 2. Get Information about Variables of Interest ---
varInfo = struct();
for ivar = 1:nVars
    currentVarName = varsOfInterest{ivar};
    
    % --- FIX: Check if the field exists before trying to access it ---
    if ~isfield(trials, currentVarName)
        error('getBinnedSpikeCounts:MissingField', ...
              'The input ''trials'' struct is missing the required field: "%s". Please check varsOfInterest.', currentVarName);
    end
    
    varInfo(ivar).name = currentVarName;
    varInfo(ivar).dimension = ivar;
    
    trial_vals = [trials.(currentVarName)];
    
    if isempty(options.uniqueVals)
        varInfo(ivar).uniqueVals = unique(trial_vals(~isnan(trial_vals)));
    else
        varInfo(ivar).uniqueVals = options.uniqueVals{ivar};
    end
    
    varInfo(ivar).nVals = numel(varInfo(ivar).uniqueVals);
    
    % Store all trial values for later lookup
    vars(ivar).trialVals = trial_vals;
end

% --- 3. Bin Spikes for Each Condition ---
allVarCombs = combvec(varInfo(:).uniqueVals);
trialVarsCat = vertcat(vars(:).trialVals);
nVarCombs = size(allVarCombs,2);
binEdgesVector = options.intervalStart:options.binSpacing:options.intervalEnd;
nBins = numel(binEdgesVector);

% Initialise cond struct array
cond(nVarCombs).varVals = [];
cond(nVarCombs).trials = [];
cond(nVarCombs).spikeCounts =[];

for ivarcomb = 1:nVarCombs
    currentVarComb = allVarCombs(:,ivarcomb)';
    cond(ivarcomb).varVals = currentVarComb;
    
    % Find trials that match the current combination of variable values
    matching_trials_idx = ismember(trialVarsCat', currentVarComb, 'rows');
    cond(ivarcomb).trials = trials(matching_trials_idx);
    
    nReps = numel(cond(ivarcomb).trials);
    
    % Initialise spike counts array with size (nBins, nReps, nUnits)
    cond(ivarcomb).spikeCounts = NaN * ones(nBins-1, nReps, nUnits); % -1 for bins vs edges
    
    if nReps > 0
        for itrial = 1:nReps
            trialStartTime = cond(ivarcomb).trials(itrial).(options.trialStartField);
            binEdges = trialStartTime + binEdgesVector;
            
            for iunit = 1:nUnits
                spike_times = units(iunit).(options.spikeTimesField);
                cond(ivarcomb).spikeCounts(:,itrial,iunit) = histcounts(spike_times, binEdges);
            end
        end
    else % if there are no trials
        cond(ivarcomb).spikeCounts = [];
    end
end

% --- 4. Add Binned Spike Field to Units Struct ---
outputShape = [varInfo.nVals];
if numel(outputShape) == 1
    outputShape = [1, outputShape];
end

for iunit = 1:nUnits
    % For each unit, create a cell array where each cell corresponds to a condition
    tempCell = cell(1, nVarCombs);
    for ivarcomb = 1:nVarCombs
        if ~isempty(cond(ivarcomb).spikeCounts)
            % Extract the slice of spike counts for the current unit
            tempCell{ivarcomb} = cond(ivarcomb).spikeCounts(:, :, iunit);
        else
            tempCell{ivarcomb} = [];
        end
    end
    
    % Reshape this cell array into the final multi-dimensional format
    units(iunit).(options.allSpikesField) = reshape(tempCell, outputShape);
end

end

