function mih = AOIfinder(varargin)
% mih = AOIfinder(analData)
%
% Spits out a mat to be loaded into flimscroll2 so we can determine the
% integraeted intensity within spot AOIs when the nuke is lacking a
% detectable spot
%
% INPUTS
%   analData == from flimscroll2
%
% OUTPUTS
%   noSpotMatrix == n x 6 mat featuring the nuke, frame, and location of
%   AOIs associated with nukes that contain spots in frames other than the
%   ones recored here
%       EX: [1.time(s) 2.frameNum 3.nucleus 4.binarySpot 5.noSpotXpos 6.noSpotYpos 7.relativeFrameNumber]
%
% USAGE
%   noSpotMatrix = AOIfinder(analData);

% Timothy Harden 20201020

analData = varargin{1};
spotsMat = analData.spotsMat;
nukesMat = analData.nukesMat;

% AllSpots = varargin{2};

% get all the nukes that have spots that need to be filled in:
nukesV = unique(spotsMat(:,3));

% go through each nuke, find frames where it doesn't have a spot, fill in a
% location where that spot may be:
n = size(nukesV,1);
frameNum = max(nukesMat(:,2)) - min(nukesMat(:,2)) + 1; % Number of frames we've integrated over. + 1 for fence post error
noSpotMat =[];
for i = 1:n
    % infor for nuke i:
    nukeNum = nukesV(i);
    singleNukeMat = nukesMat(nukesMat(:,3) == nukeNum, :); % [time frame nukeNum binary 5.nukeXpos 6.nukeYpos]
    % location of nuke i spots
    spotLocs = spotsMat(spotsMat(:,3) == nukeNum, :);
    % go through each frame and if there is no spot add one. eventually add
    % it to an analdata. then load this into flimscroll and make a button
    % to integrate the AOIs without a spot
    for j = 1:frameNum
        if singleNukeMat(j,4) == 0 %if there is no spot at this time...
            % find the position of the closest spot (in time) relative to the nuke position:
            % get the frame
            frame = singleNukeMat(j,2);
            % subtract the frame from the list of frames in which there is a spot for this nuke
            nearestV = abs(spotLocs(:,2) - frame);
            % get the index for the nearest frame with a spot. 
            [~,indx] = min(nearestV);
            % get the position of that nearest-in-time spot relative to the nucleus:
            frameNukePos = singleNukeMat(singleNukeMat(:,2) == spotLocs(indx,2),5:6); % nuke position for that frame in which there is  a spot
            relSpotPos = frameNukePos - spotLocs(indx,5:6);
            % Now add that position to the nuke position of the current frame that lacks a spot
            noSpotPos = singleNukeMat(j,5:6) - relSpotPos; 
            %put it all together
            noSpotMat = [noSpotMat; singleNukeMat(j,1:3) 0 noSpotPos singleNukeMat(j,11)]; % [1.time(s) 2.frameNum 3.nucleus 4.binarySpot 5.noSpotXpos 6.noSpotYpos 7.relativeFrameNumber]
        end
    end
end
mih = noSpotMat;
            