function [mih,plotHandle] = initialFractionBin(analData1,analData2,bins,binsOfInterest,fitInputs,tm,tx,fignum)
% [mih,plotHandle] = initialFractionBin(analData1,analData2,bins,binsOfInterest,fitInputs,tm,tx,fignum)    cia,nukesMat,normalization,fitInputs,tm,tx,fignum
% 
% As for initialFraction5/6, but to fit a subset of the events (ant/post,
% etc... as specified in t1p54h.m. Bc of details we initially have to treat each
% replicate independently
%
% INPUTS
%   bins == a vector covering the AP length of interest: eg bins = [-0.20:0.02:0.20];
%   binsOfInterest == the indicies of the bins that you want to select events from: eg boi = 5:8;
%   fitInputs == two member vector. Optional (leave empty). [ap Tau0]; note that Af = ap^2/(1+ap^2)
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
% Timothy Harden 20201106

%initialize some variables:
ciaCell{1} = analData1.cia;
ciaCell{2} = analData2.cia;
nukesMat = [analData1.nukesMat; analData2.nukesMat]; %purely to determine length of experiment


% first place the nukes into bins:
% rep 1
out1 = positionDistribution(analData1,bins); 
% rep 2
out2 = positionDistribution(analData2,bins);
% get the number of nukes in each bin:
nBinTotals = out1.totalBin + out2.totalBin;
% replace any zeros with the mean nuke num
for i = 1:length(nBinTotals)
    mn = mean(nBinTotals(nBinTotals ~= 0));
    if nBinTotals(i) == 0 
        nBinTotals(i) = round(mn);
    end
end
%determine nukes per bin:
nukesPerBin = out1.totalBin + out2.totalBin;
% replace any unwanted zeros with the mean nuke num
for i = 1:length(nukesPerBin)
    mn = mean(nukesPerBin(nukesPerBin ~= 0));
    if nukesPerBin(i) == 0 
        nukesPerBin(i) = round(mn);
    end
end
%for each rep, pull out the nukes that we want:
eventPosCell{1} = out1.eventPos; %[1.nukeNumber 2.nukePostion 3.binAssignment]
eventPosCell{2} = out2.eventPos;

comboCia = [];
for i = 1:2
    ciaMat = ciaCell{i};
    eventPosMat = eventPosCell{i};
    for j = 1:length(binsOfInterest)
        % get the nukes that belong in bin j
        logi1 = eventPosMat(:,3) == binsOfInterest(j);
        % if there are any nukes:
        if sum(logi1) ~= 0
            nukeNums = eventPosMat(logi1,1);
        % add these to a new cia, where we dont care about nuke number
            for k = 1:length(nukeNums)
                logi2 = ciaMat(:,1) == nukeNums(k);
                comboCia = [comboCia;ciaMat(logi2,:)];
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%
% cia : '1.nucleus 2.first,middle,or Last(= -3,1,or 3) 3.frameStartRelativeToNCstart 4.frameEnd 5.deltaFrames 6.deltaTime(sec) 7.startTimeRelativeToNCstart'
% nukes: '[1.time(s) 2.frameNum 3.nucleus 4.binarySpot 5.nukeXpos 6.nukeYpos 7.spotGaussAmp 8.spotSigma 9.spotOffset 10.spotIntegratedIntensity 11.relativeFrameNumber 12.NC 13.nuclearAPpos]'

%grab the first event
firstTimes = comboCia(comboCia(:,2) == -3,7); %a vector
%make the first event occur at t = 0;
first = (firstTimes - min(firstTimes) + 1); %harden 202009: you can't have your first event start at zero 
% if you want to get rid of zeros:
% first(first == 0) = 1;
% nordmalize to total number of nukes in the bins of interest:
mx = sum(nukesPerBin(binsOfInterest));
%cumulative intervals:
sints = [];
for i = 0:ceil(max(nukesMat(:,1)) - min(firstTimes)) + 1 %here we take the length of the experiment (during NC14) from the nukes mat and subtract the time of the 1st even
    logi = first(:) < i;
    sints = [sints;i sum(logi)];
end
%normailization?
sints(:,2) = sints(:,2)./mx;
%plot fit:
figure(fignum);
mih.ints = first;
mih.n = length(first);
mih.totNukes = mx;

% plot data
% plotHandle = stairs(sints(:,1) - 150,sints(:,2));shg %TSing
hold on;plotHandle = stairs(sints(:,1),sints(:,2),'LineWidth',1);shg

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
    plot(fitTs,fit,'k--');shg
end
try
    mih.Af = Af;
    mih.k0 = k0;
end
    
     
