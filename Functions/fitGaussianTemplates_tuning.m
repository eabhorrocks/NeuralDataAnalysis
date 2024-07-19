function [bestParams, char, bestR2] = fitGaussianTemplates_tuning(spikeCellArray,customOffset, plotFlag)


% do basic analysis of tuning to get initial guesses for mu at least.
[maxVal, maxidx] = max(cellfun(@mean, spikeCellArray));
[minVal, minidx] = min(cellfun(@mean, spikeCellArray));
rangeVal = range(cellfun(@mean, spikeCellArray));


%% template matching gaussian fits
gaussFun =  @(params,xdata) params(1) + params(2).*exp(-(((xdata-params(3)).^2)/(2*(params(4).^2))));
opts = optimset('Display','off');

%% check shape of spikeCellArray and the cells, and reshape if necessary

[a,b] = size(spikeCellArray);

if sum(unique([a,b])>1) > 1
    error('spikeCellArray be of size [n 1], or can be reshaped from [1 n]')
end 

if b>a % need to reshape
    spikeCellArray = spikeCellArray';
end

clear a b

[a,b] = size(spikeCellArray{1});

if sum(unique([a,b])>1) > 1
    error('spikeCellArray be of size [n 1], or can be reshaped from [1 n]')
end

if b>a % need to reshape
    spikeCellArray = cellfun(@(x) x', spikeCellArray,'UniformOutput',false);
end


%% get x and y train vals

spikes = spikeCellArray;


xtrainvals =[];
for istim = 1:numel(spikes)
    xtrainvals = [xtrainvals; repelem(istim, numel(spikes{istim}),1)];
end

ytrainvals = vertcat(spikes{:});




%% fit different template models get spike counts for this unit

% low pass (decay)
% positive amplitude
% [base, amp, mu, sigma];
lb_low1 = [-200, 0.1, 1-customOffset*10, 0.7];
ub_low1 = [200, 200, 1+customOffset 50/10];
x0_low1 = [0, rangeVal, 1-customOffset*2, 10/10];

% negative amplitude
%lb_low2 = [0, -200, 5.5, 0.7];
lb_low2 = [-200, -200, numel(spikes)-customOffset, 0.7];
ub_low2 = [200, -0.1, numel(spikes)+customOffset*10 50/10];
x0_low2 = [0, -rangeVal, numel(spikes)+customOffset*2, 10/10];

[param_out(1,:), resnorm(1)] = lsqcurvefit(gaussFun,x0_low1,xtrainvals,ytrainvals,lb_low1,ub_low1,opts);
[param_out(2,:), resnorm(2)] = lsqcurvefit(gaussFun,x0_low2,xtrainvals,ytrainvals,lb_low2,ub_low2,opts);


% highpass (rise)
% positive amplitude
% [base, amp, mu, sigma];
%lb_high1 = [0, 0.1, 5.5, 0.7];
lb_high1 = [-200, 0.1, numel(spikes)-customOffset, 0.7];
ub_high1 = [200, 200, numel(spikes)+customOffset*10, 50/10];
x0_high1 = [0, rangeVal, numel(spikes)+customOffset*2, 10/10];

% negative amplitude
lb_high2 = [-200, -200, 1-customOffset*10, 0.7];
%ub_high2 = [200, -0.1, 1.5, 10];
ub_high2 = [200, -0.1, 1+customOffset, 50/10];
x0_high2 = [0, -rangeVal, 1-customOffset*2, 10/10];

[param_out(3,:), resnorm(3)] = lsqcurvefit(gaussFun,x0_high1,xtrainvals,ytrainvals,lb_high1,ub_high1,opts);
[param_out(4,:), resnorm(4)] = lsqcurvefit(gaussFun,x0_high2,xtrainvals,ytrainvals,lb_high2,ub_high2,opts);



% bandpass (peak)
lb_band = [-200, 0.1, 1+customOffset, 0.7];
ub_band = [200, 200, numel(spikes)-customOffset, 3];
x0_band = [0, rangeVal, maxidx, 1];

[param_out(5,:), resnorm(5)] = lsqcurvefit(gaussFun,x0_band,xtrainvals,ytrainvals,lb_band,ub_band,opts);


% inverse (peak)
%lb_inv = [0, -200, 2, 0.7];
%ub_inv = [200, -0.1, 5 10];
lb_inv = [-200, -200, 1+customOffset, 0.7];
ub_inv = [200, -0.1, numel(spikes)-customOffset,3];
x0_inv = [0, -rangeVal, minidx, 1];

[param_out(6,:), resnorm(6)] = lsqcurvefit(gaussFun,x0_inv,xtrainvals,ytrainvals,lb_inv,ub_inv,opts);




idx_idx = 1;
done = false;

% get quality of fits in order
[~, fitIdx] = mink(resnorm, 7);

while ~done
    
bestIdx = fitIdx(idx_idx);
bestParams = param_out(bestIdx,:);

if bestIdx == 1 || bestIdx == 2
    char = 1; % decay
    
elseif bestIdx == 3 || bestIdx == 4
    char = 2; % rise
    
elseif bestIdx == 5 
    char = 3; % peak
    evalOfFun = feval(gaussFun,[param_out(bestIdx,:)], 1:0.1:numel(spikes));
    [~, ~, ~, prom] = findpeaks(evalOfFun);
    if prom < rangeVal/3 %0.2
        idx_idx = idx_idx+1;
        %disp('skip');
        continue
    end
    promReq = rangeVal/3; %normally 3
elseif bestIdx == 6
    char = 4; % trough
    evalOfFun = -1.*feval(gaussFun,[param_out(bestIdx,:)], 1:0.1:numel(spikes));
    %evalOfFun_norm = (evalOfFun-minVal)./(maxVal-minVal);
    [~, ~, ~, prom] = findpeaks(evalOfFun);    
    if prom < rangeVal/3 %0.2
        idx_idx = idx_idx+1;
        %disp('skip');
        continue
    end
    promReq = rangeVal/3;
end

done = true;


end

SS_tot = sum((ytrainvals - mean(ytrainvals)).^2);
bestR2 = 1-(resnorm(bestIdx)/SS_tot);



if plotFlag
    hold off
    plot(cellfun(@mean, spikes));
    hold on
    evalOfFun = feval(gaussFun,[param_out(bestIdx,:)], 1:0.1:numel(spikes));
    plot(1:0.1:numel(spikes), evalOfFun, 'r-')
    title([char, param_out(bestIdx,4)])
    pause = 1;
    
    
end