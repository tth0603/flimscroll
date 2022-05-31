function mih = binnedCiaMaker(analData1,analData2,bins,binsOfInterest)
% mih = binnedCiaMaker(analData1,analData2,bins,binsOfInterest)
% 
% will make a cia from two analData files with only events within
% binsOfInterest
%
% INPUTS
%   bins == a vector covering the AP length of interest: eg bins = [-0.20:0.02:0.20];
%   binsOfInterest == the indicies of the bins that you want to select events from: eg boi = 5:8;
%
% OUTPUS
%   MIH = a cia
%
% Timothy Harden 20201111

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
for i = 1:2 % one for each analData
    ciaMat = ciaCell{i};
    eventPosMat = eventPosCell{i};
    for j = 1:length(binsOfInterest)
        % get the nukes that belong in bin j
        logi1 = eventPosMat(:,3) == binsOfInterest(j);
        % if there are any nukes:
        if sum(logi1) ~= 0
            nukeNums = eventPosMat(logi1,1);
        % add these to a new cia, 
            for k = 1:length(nukeNums)
                logi2 = ciaMat(:,1) == nukeNums(k);
                addOn = ciaMat(logi2,:);
                if i == 2
                    addOn(:,1) = addOn(:,1) + 1000; % tack on a number so we dont confuse nukes from two different experiments
                end
                comboCia = [comboCia;addOn];
            end
        end
    end
end
mih = comboCia;