function mih = idleCiaMaker2(cia,acqnRate)
% out = idleCiaMaker(cia)
%
% this will pull out all the idle periods from a cia and spit out a cia,
% all rose of which correspond to inter-event idle periods. 
%
% v2 includes the time between the end of the last transcription event for
% a given nucleus and the end of the last trxn event across teh embryo
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

% get the time of the end of the last transcritpion event in this embryo:
endTimes = sum([cia(:,6) cia(:,7)]);
finalActiveTime = max(endTimes);
finalActiveFrame = max(cia(:,4));

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
        idleMat = zeros([eventNum,cols]); 
        for j = 1:eventNum
            idleMat(j,1) = nukeNum;
            if j == 1
                idleMat(j,2) = -2; % binary classification of the first idle period in the record
                idleMat(j,3) = mat(j,4); % the first frame of this idle period is the last frame of the corresponding trxn event
                idleMat(j,4) = mat(j + 1,3); % the last frame of this period is the first frame of the next trxn event
                idleMat(j,5) = idleMat(j,4) - idleMat(j,3); % delta frames
    %             idleMat(j,6) = idleMat(j,5)*acqnRate - acqnRate/2; % another way. this aligns with what is spit out of frequency_dweelFli7, but is not correct. 
                idleMat(j,6) = idleMat(j,5)*acqnRate; % delta time. 
                idleMat(j,7) = mat(j,6) + mat(j,7); % start time relative to the start of the NC. The start time plus the delta time of the preceeding trxn event
            elseif j < eventNum
                idleMat(j,2) = 0; % binary classification of a 'low' event 
                idleMat(j,3) = mat(j,4); % the first frame of this idle period is the last frame of the corresponding trxn event
                idleMat(j,4) = mat(j + 1,3); % the last frame of this period is the first frame of the next trxn event
                idleMat(j,5) = idleMat(j,4) - idleMat(j,3); % delta frames
    %             idleMat(j,6) = idleMat(j,5)*acqnRate - acqnRate/2; % another way. this aligns with what is spit out of frequency_dweelFli7, but is not correct. 
                idleMat(j,6) = idleMat(j,5)*acqnRate; % delta time. 
                idleMat(j,7) = mat(j,6) + mat(j,7); % start time relative to the start of the NC. The start time plus the delta time of the preceeding trxn event
            elseif j == eventNum           
                % is this the last active trxn event in this embryo?
                if mat(j,4) == finalActiveFrame
                    idleMat(j,1:7) = 0;
                else
                    idleMat(j,2) = 2; % binary classification fo the last idle period in the record, spanning the end of the last active trxn for this nucleus, and the end of the last transcription event across the embryo
                    idleMat(j,3) = mat(j,4); % the first frame of this idle period is the last frame of the corresponding trxn event
                    idleMat(j,4) = finalActiveFrame; % the last frame of this period is the first frame of the next trxn event
                    idleMat(j,5) = idleMat(j,4) - idleMat(j,3); % delta frames
                    idleMat(j,6) = idleMat(j,5)*acqnRate; % delta time. 
                    idleMat(j,7) = mat(j,6) + mat(j,7); % start time relative to the start of the NC. The start time plus the delta time of the preceeding trxn event
                end
            end
        end
    elseif eventNum == 1
        idleMat = zeros([1,cols]);
        % is this the last active trxn event in this embryo?
        if mat(1,4) == finalActiveFrame
            idleMat(1,1:7) = 0;
        else
            idleMat(1,1) = nukeNum;
            idleMat(1,2) = 2; % binary classification fo the last idle period in the record, spanning the end of the last active trxn for this nucleus, and the end of the last transcription event across the embryo
            idleMat(1,3) = mat(1,4); % the first frame of this idle period is the last frame of the corresponding trxn event
            idleMat(1,4) = finalActiveFrame; % the last frame of this period is the first frame of the next trxn event
            idleMat(1,5) = idleMat(1,4) - idleMat(1,3); % delta frames
            idleMat(1,6) = idleMat(1,5)*acqnRate; % delta time. 
            idleMat(1,7) = mat(1,6) + mat(1,7); % start time relative to the start of the NC. The start time plus the delta time of the preceeding trxn event
        end
    end
    idleCia = [idleCia; idleMat];
end
% get rid of the row pertaining to the (non existant) idle period after the
% end of teh last trxn event in the embryo:
idleCia(idleCia(:,1) == 0,:) = [];
mih = idleCia;
    