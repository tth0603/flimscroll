function mih = finalIdle(varargin)
%
% This script reports a vector of the intervals between the end of the last
% active trxn event within each nuclei and the end of last trxn event
% observed across the entire embryo. A shitty proxy for the start of
% gastrulation. 
%
% see also finalIdle2 for a better, albeit more time intensive, way. 
%
% INPUTS:
% varargin == a bunch of analDatas. Not center cias. Not dropout cias
%
% OUTPUTS:
% a vector of times you can plot with survivalPlotFromVector. These are the
% times between the last active transcription event and the "start" of
% gastrulation.

% initialize:
repNum = size(varargin,2);

% now define the center of the stripe (this could be an input, but we've
% never changed this so I just wrote it in:
bins = [-0.20:0.02:0.20]; % define the bins
boiC = 10:11; % for the center bins

% Now get the intervals:
endEvents = [];
for i = 1:repNum
    % grab the raw cia for this rep:
    cia = varargin{i}.cia;

%     % too use the last frame as the end of the experimetn, not the end of
%     % the last active trxn interval:
%     nukesMat = varargin{i}.nukesMat;
%     finalFrame = max(nukesMat(:,11)); % this is the last frame relative to the start of NC14
%     logi = nukesMat(:,3) == 1;
%     timeV = nukesMat(logi,1);
%     meanFrameInt = mean(diff(timeV(:,1))); % this is the average time between frames (in seconds)
%     finalTime = finalFrame*meanFrameInt;
% 
% the above is gobbilty gook, or at leasrt it won't work. We need a better
% measure of the start of gastrulation than the end of the experiment. 

    %find the time that the last active trxn time ends, use this as the end
    %of the experiment:
    lastEvents = cia(cia(:,2) == 3,:);
    lastTimes = lastEvents(:,7) + lastEvents(:,6); % time they begin rel to NC 14 start plus the event length
    lastTime = max(lastTimes);
    % now look at only the center nukes:
    centerCia = binnedCiaMaker3(varargin{i},bins,boiC);
    % find the number of nukes in the center:
    nukes = unique(centerCia(:,1));
    len = length(nukes);
    % now find the end idle period for each nuke in the center:
    for j = 1:len
        mat2 = centerCia(centerCia(:,1) == nukes(j),:);
        endInt = lastTime - (mat2(end,7) + mat2(end,6));
%         endInt = finalTime - (mat2(end,7) + mat2(end,6)); % this is looking at the end of the experiment interval times. It's garbage
        if endInt == 0
            endInt = [];
        end
        endEvents = [endEvents; endInt];
    end
end

mih = endEvents;