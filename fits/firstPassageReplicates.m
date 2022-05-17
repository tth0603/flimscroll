function [v1,v2] = firstPassageReplicates(analData1,analData2,bins,binsOfInterest)
% [mih,plotHandle] = initialFractionBin(analData1,analData2,bins,binsOfInterest,fitInputs,tm,tx,fignum)    cia,nukesMat,normalization,fitInputs,tm,tx,fignum
% 
% As for initialFractionBin2, spits out vectors of events for comparison
% wiht KS test
%
% INPUTS
%   bins == a vector covering the AP length of interest: eg bins = [-0.20:0.02:0.20];
%   binsOfInterest == the indicies of the bins that you want to select events from: eg boi = 5:8;
%
% OUTPUS
%   out1, out2 = each of vector of first passage times from the center
%   nukes of each replicate
%
% [v1, v2] = firstPassageReplicates(analData1,analData2,bins,binsOfInterest)
%cia descrip: '1.nucleus 2.first,middle,or Last(= -3,1,or 3) 3.frameStart 4.frameEnd 5.deltaFrames 6.deltaTime(sec) 7.startTime';
%
% Timothy Harden 20210616

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

% for rep 1:
mat1 = [];
ciaMat = ciaCell{1};
eventPosMat = eventPosCell{1};
for j = 1:length(binsOfInterest)
    % get the nukes that belong in bin j
    logi1 = eventPosMat(:,3) == binsOfInterest(j);
    % if there are any nukes:
    if sum(logi1) ~= 0
        nukeNums = eventPosMat(logi1,1);
    % add these to a new cia, where we dont care about nuke number
        for k = 1:length(nukeNums)
            logi2 = ciaMat(:,1) == nukeNums(k);
            mat1 = [mat1;ciaMat(logi2,:)];
        end
    end
end
% output a vector:
% v1 = mat1(mat1(:,2) == -3,7); %a vector, uses time
v1 = mat1(mat1(:,2) == -3,3); %a vector, uses frames

% for rep 1:
mat2 = [];
ciaMat = ciaCell{2};
eventPosMat = eventPosCell{2};
for j = 1:length(binsOfInterest)
    % get the nukes that belong in bin j
    logi1 = eventPosMat(:,3) == binsOfInterest(j);
    % if there are any nukes:
    if sum(logi1) ~= 0
        nukeNums = eventPosMat(logi1,1);
    % add these to a new cia, where we dont care about nuke number
        for k = 1:length(nukeNums)
            logi2 = ciaMat(:,1) == nukeNums(k);
            mat2 = [mat2;ciaMat(logi2,:)];
        end
    end
end
% output a vector:
% v2 = mat2(mat2(:,2) == -3,7); %a vector, use time
v2 = mat2(mat2(:,2) == -3,3); %a vector, use frames
