function [reliability, tuningVar, allVar] = calcReliability(spikeCountCellArray)
            
% from Christensen and Pillow, 2022 N. Comms
tuningVar = var(cellfun(@mean, spikeCountCellArray));
allVar = var(cat(1,spikeCountCellArray{:}));

reliability = tuningVar/allVar;
end

