function [mib,plotHandle] = initialFraction(cia,normalization,fitInputs,fignum)
% [mib,plotHandle] = initialFraction(cia,normalization,fitInputs,fignum)
%cia descrip: '1.nucleus 2.first,middle,or Last(= -3,1,or 3) 3.frameStart 4.frameEnd 5.deltaFrames 6.deltaTime(sec) 7.startTime';

%grab the first event (its length)
first = cia(cia(:,2) == -3,7); %a vector

%determine normalization
if ~isempty(normalization)
    mx = normalization;
else 
%     mx = length(first); %there is a debate in my head as to which is more
%     useful
    mx = 1;
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
    fitParams = fminsearch('twoStepFirstPassage',fitInputs,[],first);
    fitTs = 0:ceil(max(first));
    k1 = fitParams(1);
    k2 = fitParams(2);
%     k1 = fitInputs(1);
%     k2 = fitInputs(2);
    fit = 1-(k1*k2/(k1 - k2).*((1/k2).*exp(-k2.*fitTs)-(1/k1).*exp(-k1.*fitTs)));
    hold on;
    plot(fitTs,fit,'k-');shg
end
try
    mib.k1 = k1;
    mib.k2 = k2;
end
     
