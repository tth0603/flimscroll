function mih = binnedCiaMaker2(analData,bins,binsOfInterest)
% mih = binnedCiaMaker(analData1,analData2,bins,binsOfInterest)
% 
% will make a cia with only events within binsOfInterest
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
cia = analData.cia;

% first place the nukes into bins:
% rep 1
out = positionDistribution(analData,bins); 
% get the number of nukes in each bin:
nBinTotals = out.totalBin;
% replace any zeros with the mean nuke num
for i = 1:length(nBinTotals)
    mn = mean(nBinTotals(nBinTotals ~= 0));
    if nBinTotals(i) == 0 
        nBinTotals(i) = round(mn);
    end
end
%determine nukes per bin:
nukesPerBin = out.totalBin;
% replace any unwanted zeros with the mean nuke num
for i = 1:length(nukesPerBin)
    mn = mean(nukesPerBin(nukesPerBin ~= 0));
    if nukesPerBin(i) == 0 
        nukesPerBin(i) = round(mn);
    end
end
%for each rep, pull out the nukes that we want:
eventPos = out.eventPos; %[1.nukeNumber 2.nukePostion 3.binAssignment]

newCia = [];
for j = 1:length(binsOfInterest)
    % get the nukes that belong in bin j
    logi1 = eventPos(:,3) == binsOfInterest(j);
    % if there are any nukes:
    if sum(logi1) ~= 0
        nukeNums = eventPos(logi1,1);
    % add these to a new cia, 
        for k = 1:length(nukeNums)
            logi2 = cia(:,1) == nukeNums(k);
            addOn = cia(logi2,:);
            newCia = [newCia;addOn];
        end
    end
end

mih = newCia;