function mih = positionDistribution(analData,bins)
% out = positionDistribution(analData,bins)
%
% This will find the center of the stripe pattern in the embryo
% corresponding to analData, then bin the nukes around that center. This is
% because our method of absolute AP position measurement isn't as precise as we'd
% like. 
%
% INPUTS
%   analData is the data strucure
%   bins == a vector covering the AP length of interest with the center of the stripe pattern at 0: eg bins = [-0.20:0.02:0.20];
%
% OUTPUTS
%   out.dist == distribution of fraction of nukes within each bin from
%   'bins' that ever have a spot
%   out.eventPos == two column vector of [1.nukeNumber 2.nukePostion 3.binAssignment] for
%   nukes with an event
%   out.noEventNukePos == as above, but nukes without an event.
%
% Harden 2020

% initiate variables:
nukeMat = analData.nukesMat;
ciaMat = analData.cia;
% for nukes with an event:
eventNukes = nukeMat(~isnan(nukeMat(:,7)),:); %grab only rows with a detected spot
nukeNums = unique(eventNukes(:,3)); %the numer of unique nuclei with an event
% number of nukes:
n = size(nukeNums,1);
%now go back and get the info for these nukes for all time points
nukesWithEvent = [];
for i = 1:n
    logi = nukeMat(:,3) == nukeNums(i);
    nukesWithEvent = [nukesWithEvent; nukeMat(logi,:)];
end
% find the mid point:
timepoints = max(nukeMat(:,2)) - min(nukeMat(:,2)) + 1; %number of time points
%here is the position (mean) of each nuke:
nukePosMat = zeros(timepoints,n);
for i = 1:n
    logi = nukesWithEvent(:,3) == nukeNums(i);
    nukePosMat(:,i) = nukesWithEvent(logi,13); %all nuke positions at all time points
end
nukePosV = [nukeNums mean(nukePosMat,1)']; %the mean position of each nuke across nc 14
stripeCenter = mean(nukePosV(:,2));% the center of the stripe in terms of absolute AP position
% center the positions
nukePosV(:,2) = nukePosV(:,2) - stripeCenter;

% for nukes without an event:
allNukeNums = unique(nukeMat(:,3));
noEventNukeNums = setxor(nukeNums,allNukeNums);
% number of no event nukes:
noEventN = size(noEventNukeNums,1);
%now go back and get the info for these nukes for all time points
nukesWithOutEvent = [];
for i = 1:noEventN
    logi = nukeMat(:,3) == noEventNukeNums(i);
    nukesWithOutEvent = [nukesWithOutEvent; nukeMat(logi,:)];
end
% find the mid point:
%here is the position (mean) of each no event nuke:
noEventNukePosMat = zeros(timepoints,noEventN);
for i = 1:noEventN
    logi = nukesWithOutEvent(:,3) == noEventNukeNums(i);
    noEventNukePosMat(:,i) = nukesWithOutEvent(logi,13); %all nuke positions at all time points
end
noEventNukePosV = [noEventNukeNums mean(noEventNukePosMat,1)']; %the mean position of each nuke across nc 14
% center the positions
noEventNukePosV(:,2) = noEventNukePosV(:,2) - stripeCenter;
% now come up with some distributions:
[eN,eventBin] = histc(nukePosV(:,2),bins);
[neN,noEventBin] = histc(noEventNukePosV(:,2),bins);
% output
mih.eventBin = eN;
mih.totalBin = eN + neN;
mih.eventPos = [nukePosV eventBin]; %[1.nukeNumber 2.nukePostion 3.binAssignment]
mih.noEventNukePos = [noEventNukePosV noEventBin]; %[1.nukeNumber 2.nukePostion 3.binAssignment]

