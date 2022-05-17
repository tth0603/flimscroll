function mih = edgeFinder(nukesMat,cia,edges,split)
% out = edgeFinder(nukesMat,cia,edges,split)
%
% Takes a single cia and uses the AP position (col. 13) in nukesMat to
% split the events into two subset cias spatially: either those at the
% edges/middle of the stripe or those in the anterior/posterior
%
% INPUTS:
%   nukesMat == from analData.nukesMat
%   cia == from analData.cia
%   edges == binary. 1 to find nukes on either edge of the stripe. 0 to
%       split the stripe in two parts
%   split == 0 to 1. if EDGES this specifies the fraction of nukes we look
%       at on both edges. If else this this specifies the fraction of nukes
%       in the ANT bin. 
%
% OUTPUTS:
%   out == a structure containing two matrici members:
%   out.cia1 == for all events on the edges OR in the ANT side of the stripe
%   out.cia2 == all events in the middle OR the POST side
%   out.nukePosition
%
% USAGE:
%   out = edgeFinder(nukesMat,cia,1,0.2)
%
% Timothy Harden 20200803

nukeMat = nukesMat;
ciaMat = cia;
% get nukes with an event:
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
nukePosV = mean(nukePosMat,1); %the mean position of each nuke across nc 14
stripeCenter = mean(nukePosV); % the center of the stripe in terms of AP position
% figure(2);histogram(nukePosV,20);shg % to TS the distribution of events across the AP
% for binning nukes:
sNukePosV = sort(nukePosV); %order nuke positions from least to greatest
if edges == 1
    antEdge = sNukePosV(round(split*n)); %get the index that pertains to the SPLIT-th percentile, then pluck that out of the list of nuke positions
    postEdge = sNukePosV(round((1 - split)*n)); %get the index that pertains to the (1 - SPLIT)-th percentile, then pluck that out of the list of nuke positions
    %the list of nukes in the edges:
    nukeList1 = []; % these are the edge nukes
    nukeList2 = []; % these are the mid nukes
    for i = 1:n
        if nukePosV(i) < antEdge || nukePosV(i) > postEdge % for edge nukes they must be less than the antEdge AP value OR greater than the postEdge AP value
            nukeList1 = [nukeList1; nukeNums(i)];
        else 
            nukeList2 = [nukeList2; nukeNums(i)];
        end
    end    
else 
    antBin = sNukePosV(round(split*n));
    %the list of nukes on each side of the stripe:
    nukeList1 = []; % these are the ant nukes
    nukeList2 = []; % these are the post nukes
    for i = 1:n
        if nukePosV(i) < antBin % for ANT nukes they must be less than the antBin AP value 
            nukeList1 = [nukeList1; nukeNums(i)];
        else 
            nukeList2 = [nukeList2; nukeNums(i)];
        end
    end    
end
% now pull out the corresponding events from the cia
L1 = size(nukeList1,1);
L2 = size(nukeList2,1);
cia1 = [];
cia2 = [];
for i = 1:L1
    logi1 = ciaMat(:,1) == nukeList1(i);
    cia1 = [cia1; ciaMat(logi1,:)];
end
for i = 1:L2
    logi2 = ciaMat(:,1) == nukeList2(i);
    cia2 = [cia2; ciaMat(logi2,:)];
end 
mih.cia1 = cia1;
mih.cia2 = cia2;