function [mih,plotHandle] = initialFractionBin2(varargin)
% [mih,plotHandle] = initialFractionBin(analData1,analData2,...,bins,binsOfInterest,fitInputs,tm,tx,fignum)    
% 
% As for initialFractionBin, but the curves are not x-axis shifted. 
% 
% INPUTS
%   bins == a vector covering the AP length of interest: eg bins = [-0.20:0.02:0.20];
%   binsOfInterest == the indicies of the bins that you want to select events from: eg boi = 5:8;
%   fitInputs == two member vector. Optional (leave empty). [ap Tau0]; note that Af = ap^2/(1+ap^2)
%   tm (tx) == minimum (maximum) time resolution. For fits. Optional (leave empty)
%
% OUTPUS
%   MIH = a four member struc with members: Af (scaler), k0 (scaler), n
%   (scaler, number of events), and ints (a vector of events after they
%   have been shifted so first event is at the origin).
%   plotHandle = to use with SET e.g. set(plotHandle,'Color',colr{i},'LineWidth',1)
%
%cia descrip: '1.nucleus 2.first,middle,or Last(= -3,1,or 3) 3.frameStart 4.frameEnd 5.deltaFrames 6.deltaTime(sec) 7.startTime';
%
% Timothy Harden 20201106

% number of replicates:
numberOfReps = length(varargin) - 6;
nukesMat = []; %purely to determine length of experiment
for i = 1:numberOfReps
    analData{i} = varargin{i};
    ciaCell{i} = analData{i}.cia;
    nukesMat = [nukesMat;analData{i}.nukesMat];
end

% initialize other variables:
bins = varargin{numberOfReps + 1};
binsOfInterest = varargin{numberOfReps + 2};
fitInputs = varargin{numberOfReps + 3};
tm = varargin{numberOfReps + 4};
tx = varargin{numberOfReps + 5};
fignum = varargin{numberOfReps + 6};

% place the nukes into bins:
for i = 1:numberOfReps
    out{i} = positionDistribution(analData{i},bins);
end
% get teh number of nukes in each bin:
nukesPerBin = zeros(size(out{1}.totalBin,1),1);
for i = 1:numberOfReps
    nukesPerBin = nukesPerBin + out{i}.totalBin;
end
% replace any zeros with the mean nuke num:
for i = 1:size(nukesPerBin,1)
    mn = mean(nukesPerBin(nukesPerBin ~= 0));
    if nukesPerBin(i) == 0
        nukesPerBin(i) = round(mn);
    end
end
% for each rep, pull out the nukes we want:
for i = 1:numberOfReps
    eventPosCell{i} = out{i}.eventPos; %[1.nukeNumber 2.nukePostion 3.binAssignment]
end

comboCia = [];
for i = 1:numberOfReps
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
first = comboCia(comboCia(:,2) == -3,7); %a vector
% nordmalize to total number of nukes in the bins of interest:
mx = sum(nukesPerBin(binsOfInterest));

%cumulative intervals:
sints = [];
for i = 0:ceil(max(nukesMat(:,1))) + 1 %here we take the length of the experiment (in seconds; during NC14; + 1 in case events last until the end of acquisition) from the nukes mat
    logi = first(:) < i;
    sints = [sints;i sum(logi)];
end
%normailization:
sints(:,2) = sints(:,2)./mx;
%plot fit:
figure(fignum);
mih.ints = first;
mih.n = length(first);
mih.totNukes = mx; % the total number of nukes within the bins of interest
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
    plot(fitTs,fit,'k-','LineWidth',1);shg
end
% plot data:
% plotHandle = stairs(sints(:,1) - 150,sints(:,2));shg %TSing
hold on;plotHandle = stairs(sints(:,1),sints(:,2),'LineWidth',1);shg
% if you used a fit, report the params:
try
    mih.Af = Af;
    mih.k0 = k0;
end
    
     
