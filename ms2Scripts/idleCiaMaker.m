function mih = idleCiaMaker(cia,acqnRate)
% out = idleCiaMaker(cia)
%
% this will pull out all the idle periods from a cia and spit out a cia,
% all rose of which correspond to inter-event idle periods. This ignores
% periods from the end of the last trxn event until gastrulation. 
%
% INPUTS:
%   cia == cia
%   acqnRate == the freq of acquisition. A scaler. 
%       EX: 2 frames per minute is '30' for one frame every 30s. 
%
% Harden 20220404

% pull out the nukes: 
nukes = unique(cia(:,1));
len = length(nukes);

idleCia = [];
% pull out idle period info:
for i = 1:len
    % get teh interval infor for nuke i:
    nukeNum = nukes(i);
    logi = cia(:,1) == nukeNum;
    mat = cia(logi,:);
    % if the nuke has multiple events, pull out the idle period(s)
    [eventNum,cols] = size(mat);
    if eventNum > 1
        idleMat = zeros([eventNum - 1,cols]); % -1 bc there is one less idle period than event
        for j = 1:eventNum - 1
            idleMat(j,1) = nukeNum;
            if j == 1
                idleMat(j,2) = -2; % binary classification of the first idle period in the record
            else
                idleMat(j,2) = 0; % binary classification of a 'low' event 
            end
            idleMat(j,3) = mat(j,4); % the first frame of this idle period is the last frame of the corresponding trxn event
            idleMat(j,4) = mat(j + 1,3); % the last frame of this period is the first frame of the next trxn event
            idleMat(j,5) = idleMat(j,4) - idleMat(j,3); % delta frames
%             idleMat(j,6) = idleMat(j,5)*acqnRate - acqnRate/2; % another way. this aligns with what is spit out of frequency_dweelFli7, but is not correct. 
            idleMat(j,6) = idleMat(j,5)*acqnRate; % delta time. 
            idleMat(j,7) = mat(j,6) + mat(j,7); % start time relative to the start of the NC. The start time plus the delta time of the preceeding trxn event
        end
        idleCia = [idleCia; idleMat];
    end
end
mih = idleCia;
    