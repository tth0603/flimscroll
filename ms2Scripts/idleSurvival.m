function mih = idleSurvival(cia)
% mih = idleSurvival(cia)
% 
% spits out a vector of idle times. This version
% excludes the times after the last event in each nuke
%
% Timothy Harden 20200904

% mat = dropout(cia);
mat = cia;
nukes = unique(mat(:,1));
len = length(nukes);
interEvents = [];
for j = 1:len
    mat2 = mat(mat(:,1) == nukes(j),:);
    intervals = diff(mat2(:,7)) - mat2(1:end-1,6) - 60; % start time of next event minus start time of this event plus this even't duration. The -60 is to get rid of the lag at the start of the survival curve
    interEvents = [interEvents; intervals];
end
mih = interEvents;