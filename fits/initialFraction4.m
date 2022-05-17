function [mih,plotHandle] = initialFraction4(cia,normalization,fitInputs,fignum)
%
% garbage. see v5
%
% [mih,plotHandle] = initialFraction4(cia,normalization,fitInputs,fignum)
% 
% shuves the first event to t = 0. fits curves to a two step fn with two
% equal steps: f(t) = Af*(k0)^2*t*exp(-t*k0). uses fn twoStepBindFitV2.m
%
% INPUTS
%   normalization == a scaler. the number we are normalizing the initial
%   bind curve to
%   fitInputs == model initial guesses. [Af k0] 
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
% Timothy Harden 2020

%grab the first event
first = cia(cia(:,2) == -3,7); %a vector
%make the first event occur at t = 0;
first = (first - min(first) + 1); %harden 202009: you can't have your first event start at zero 
% if you want to get rid of zeros:
% first(first == 0) = 1;
%determine normalization
if ~isempty(normalization)
    mx = normalization;
else 
%     mx = length(first);
    mx = 1;
end
%cumulative intervals:
sints = [];
for i = 0:ceil(max(cia(:,7))) 
% for i = 0:ceil(max(first)) 
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
tic;
% options = optimset('MaxIter',10^7,'MaxFunEvals',10^7);
if ~isempty(fitInputs)
%     fitParams = fminsearch('twoStepBindFitV2',fitInputs,options,first);
    fitParams = fminsearch('twoStepBindFitV2',fitInputs,[],first,30,3000);
    fitTs = 0:ceil(max(cia(:,7)));
    Af = fitParams(1);
    k0 = fitParams(2);
    fit = cumsum(Af*(1*k0)^2*fitTs.*exp(-fitTs*k0));
    hold on;
    plot(fitTs,fit,'k-');shg
end
toc;
try
    mih.Af = Af;
    mih.k0 = k0;
end
    
     
