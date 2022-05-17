function [mib,plotHandle] = initialFraction3(cia,normalization,fitInputs,fignum)
%
%single expl fit
%
% [mib,plotHandle] = initialFraction(cia,normalization,fitInputs,fignum)
%cia descrip: '1.nucleus 2.first,middle,or Last(= -3,1,or 3) 3.frameStart 4.frameEnd 5.deltaFrames 6.deltaTime(sec) 7.startTime';

%grab the first event
first = cia(cia(:,2) == -3,7); %a vector
%make the first event occur at t = 0;
first = first - min(first);
%determine normalization
if ~isempty(normalization)
    mx = normalization;
else 
    mx = length(first);
%     mx = 1;
end
%cumulative intervals:
sints = [];
for i = 0:ceil(max(cia(:,7)))
    logi = first(:) < i;
    sints = [sints;i sum(logi)];
end
%normailization?
sints(:,2) = sints(:,2)./mx;
%plot it:
figure(fignum);
plotHandle = stairs(sints(:,1),sints(:,2));shg
mib.ints = sints;
%add a fit:
if ~isempty(fitInputs)
    fitParams = fminsearch('expfallone_mxl',fitInputs,[],first,1,3600);
    fitTs = 0:ceil(max(first));
    k1 = fitParams;
    fit = k1*exp(-k1*fitTs);
    hold on;
    plot(fitTs,fit,'k-');shg
end
mib.k1 = k1;

    
     
