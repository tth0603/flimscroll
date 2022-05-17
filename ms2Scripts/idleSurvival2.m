function mih = idleSurvival2(cia)
% mih = idleSurvival2(cia)
% 
% spits out a vector of idle times. omits dropout times. This versiion
% includes a version of the time after the last event in a nuke. It looks
% for the time that the last active trxn interval ends, then uses that as a
% proxy for the end og gastrulation. see idleSurvival3 for a script that
% looks at the time when we turn of the camera as teh beginning of
% gastrulation. 
%
% Timothy Harden 20200904

% cia = dropout(cia);
% get the time that the last event ends (a bit trickier)\:
lastEvents = cia(cia(:,2) == 3,:);
lastTimes = lastEvents(:,7) + lastEvents(:,6); % time they begin rel to NC 14 start plus the event length
lastTime = max(lastTimes);
% nukes 
nukes = unique(cia(:,1));
len = length(nukes);
interEvents = [];
endEvents = [];
for j = 1:len
    mat2 = cia(cia(:,1) == nukes(j),:);
    %the times btwn events:
%     intervals = diff(mat2(:,7)) - mat2(1:end-1,6) - 60; % start time of next event minus start time of this event plus this even't duration. The -60 is to get rid of the lag at the start of the survival curve
    intervals = diff(mat2(:,7)) - mat2(1:end-1,6);
    %the times after events
    endInt = lastTime - (mat2(end,7) + mat2(end,6));
    if endInt == 0
        endInt = [];
    end
    interEvents = [interEvents; intervals];
    endEvents = [endEvents; endInt];
end
mih.interEvents = interEvents;
mih.endEvents = endEvents;