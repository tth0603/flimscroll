function mih = totalActiveTime(varargin)
% mih = totalActiveTime(analData1,...,analDataN)
%
% This script will spit out a vector of times. These times correspond to
% the interval between the start of the first active trxn event and the end
% of the last active trxn event within the same nuke. 
%
% It considers only nukes within the center of the stripe. 
% 
% INPUTS:
% varargin == a bunch of analDatas. Not center cias. Not dropout cias

% initialize:
repNum = size(varargin,2);

% now define the center of the stripe (this could be an input, but we've
% never changed this so I just wrote it in:
bins = [-0.20:0.02:0.20]; % define the bins
boiC = 10:11; % for the center bins

% Now get the intervals:
ints = [];
for i = 1:repNum
    % the cia for the center nukes of this rep:
    centerCia = binnedCiaMaker3(varargin{i},bins,boiC);
    
    % treat each reps individually:
    nukeNums = unique(centerCia(:,1));
    numberOfNukes = size(nukeNums,1);
    for j = 1:numberOfNukes
        % find the events that take place in nuke j
        logi = centerCia(:,1) == nukeNums(j);
        thisCia = centerCia(logi,:);
        numOfEvents = size(thisCia,1);
        % if there is only one event in this nuke:
        if numOfEvents == 1
            ints = [ints;thisCia(6)];
        else
            % get the time the first event starts:
            startTime = thisCia(1,7);
            % get teh time the last event ends:
            endTime = thisCia(end,7) + thisCia(end,6); % this is the time the last event starts plus its duration. 
            ints = [ints;endTime - startTime];
        end
    end
end
mih = ints;

        