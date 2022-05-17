function [mih,plotHandle,plotHandle2] = initialFraction7(cia,nukesMat,normalization,fitInputs,tm,tx,fignum)
% [mih,plotHandle] = initialFraction5(cia,nukesMat,normalization,fitInputs,tm,tx,fignum)
% 
% exacly like initialFraction5, but with some tweaks to plotting for a
% final Supp Fig.
%
% INPUTS
%   normalization == a scaler. the number we are normalizing the initial
%   bind curve to. if '1' we normalize to the number of events in cia, if
%   empty '[]' we do not normalize. 
%   fitInputs == Optional (leave empty). model initial guesses. [Af k0] 
%   tm == minimum time resolution
%   tx == maximum time an event can take place
%   fignum == a scaler. where the plot is placed
%
% OUTPUS
%   MIH = a four member struc with members: Af (scaler), k0 (scaler), n
%   (scaler, number of events), and ints (a vector of events after they
%   have been shifted so first event is at the origin).
%   plotHandle = to use with SET e.g. set(plotHandle,'Color',colr{i},'LineWidth',1)
%
% [mib,plotHandle] = initialFraction(cia,normalization,fitInputs,fignum)
%cia descrip: '1.nucleus 2.first,middle,or Last(= -3,1,or 3) 3.frameStart 4.frameEnd 5.deltaFrames 6.deltaTime(sec) 7.startTime';
%
% Timothy Harden 20200904

%grab the first event
firstTimes = cia(cia(:,2) == -3,7); %a vector
%make the first event occur at t = 0;
first = (firstTimes - min(firstTimes) + 1); %harden 202009: you can't have your first event start at zero 
% first = firstTimes;
% if you want to get rid of zeros:
% first(first == 0) = 1;
%determine normalization
if normalization == 1
    mx = size(first,1);
elseif ~isempty(normalization)
    mx = normalization;
else 
%     mx = length(first);
    mx = 1;
end
%cumulative intervals:
sints = [];
for i = 0:ceil(max(nukesMat(:,1)) - min(firstTimes)) + 1 %here we take the length of the experiment (during NC14) from the nukes mat and subtract the time of the 1st even
    logi = first(:) < i;
    sints = [sints;i sum(logi)];
end
%normailization?
sints(:,2) = sints(:,2)./mx;
%plot it:
figure(fignum);
% plotHandle = stairs(sints(:,1) - 150,sints(:,2));shg %TSing
plotHandle = stairs(sints(:,1),sints(:,2));shg
mih.ints = first;
mih.n = length(first);
%add a fit:
% options = optimset('MaxIter',10^7,'MaxFunEvals',10^7);
if ~isempty(fitInputs)
%     fitParams = fminsearch('twoStepBindFitV3',fitInputs,options,first);
    fitParams = fminsearch('twoStepBindFitV3',fitInputs,[],first,tm,tx,mx);
    % to change the x axis of the plot:
%     fitTs = 0:ceil(max(cia(:,7)));
    fitTs = 0:4000; 
    ap = fitParams(1);
    Af = ap^2/(1 + ap^2);
    k0 = fitParams(2);
    fit = cumsum(Af*(1/k0)^2*fitTs.*exp(-fitTs/k0));
    hold on;
    plotHandle2 = plot(fitTs,fit,'--');shg
end
try
    mih.Af = Af;
    mih.k0 = k0;
end
    
     
