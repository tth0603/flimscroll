function mib = makeInts(nukesMat)
% mib = makeInts(nukesMat)
%
% This is for use with analData.nukesMat straight outta flimscroll. It will
% give a rastergram/heatmap output of single nucleus traces.
%
% INPUTS
%   nukesMat == a nuke mat from the analData struc from flimscroll.
%       format: %[1.time(s) 2.frameNum 3.nucleus 4.binarySpot 5.nukeXpos 6.nukeYpos 7.spotGaussAmp 8.spotSigma 9.spotOffset 10.spotIntegratedIntensity 11.relativeFrameNumber 12.NC]
%
% OUPUTS
%   cumulativeIntervalsArray struc: a two member structure with cia (a
%   matrix) and its description (a string).
%   cumulativeIntervalsArrayDescription: '1.nucleus 2.first,middle,or Last(= -3,1,or 3) 3.frameStart 4.frameEnd 5.deltaFrames 6.deltaTime(sec) 7.startTime';
%
% USAGE
%   ciaStruc = makeInts(analData.nukesMat);
%
% This is now called from within flimscroll via matchMaker.
%
% Harden 2018

%get acqn frequency:
oneNuk = nukesMat(nukesMat(:,3) == 1,:); %pick all entris for first nuke
%get teh mean time between frames for this:
acqnT = mean(diff(oneNuk(:,1)));
% get the NaN's outta here 
nukesMat(isnan(nukesMat)) = 0;
% get all instances with a spot:
spotMat = nukesMat(nukesMat(:,4) > 0,:);
% get those nukes:
nukes = unique(spotMat(:,3));
%get the nukes with spots:
ints = [];
for i = 1:length(nukes)
    %get all rose for this nuke
    mat2 = nukesMat(nukesMat(:,3) == nukes(i),:);
    % a vector of binaries
    binV = mat2(:,4)';
    %find the runs of ones, where they start and stop
    isOne = [false, (binV == 1),false];
    indx = [strfind(isOne, [false, true]); ...
            strfind(isOne, [true, false]) - 1]; % [strtFrmIndx;finalFrmIndx];
    %num of events:
    noe = size(indx,2);
    for j = 1:noe
        %event type (we differentiate between first, last and final event):
        if j == 1
            eventType = -3;
        elseif j == noe(end)
            eventType = 3;
        else
            eventType = 1;
        end
        ints = [ints;nukes(i) eventType mat2(indx(1,j),2) mat2(indx(2,j),2) mat2(indx(2,j),2)-mat2(indx(1,j),2)+1 mat2(indx(2,j),1)-mat2(indx(1,j),1)+0.5*acqnT mat2(indx(1,j),1)];
        %      [   nucleusNum (-3,1,or3)    strtFrm         finalFrm              finalFrm-strtFrm+1                    finalT-strtT+0.5*acqnTime                  strtTime  ]
    end
end

mib.cia = ints;
mib.cumulativeIntervalsArrayDescription = '1.nucleus 2.first,middle,or Last(= -3,1,or 3) 3.frameStartRelativeToNCstart 4.frameEnd 5.deltaFrames 6.deltaTime(sec) 7.startTimeRelativeToNCstart';
        
