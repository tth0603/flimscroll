function mih = spotPairing(AllSpots,DL,thdx,thdy)
%
% This will associate spots to one another in time without nuclear
% segementation. Need the AllSpots struc from spot tracking and DL struc
% from making a driftlist, both in flimscroll2
%
% Make sure the drift list and the spot picking cover the same frame
% range
%
% INPUTS:
%   AllSpots == from Find all spots in flimscroll2
%   DL == from Create drift list in flimscroll 2
%   thdx/y == a threshold value for how close spots need to be. In pixels.
%       get an idea of this value via: mean(abs(DL.diffdriftlist(:,2))) [for x] and
%       mean(abs(DL.diffdriftlist(:,3))) [for y], or max(abs(DL.diffdriftlist(:,2)))...
%
% OUTPUTS:
%   A spotsMat [1.time(s) 2.frameNum 3.nucleus 4.distanceFromNucleusCntr 5.xPos 6.yPos 7.gaussAmp 8.sigma 9.offset 10.integratedIntensity 11.relativeFrameNumber 12.NC 13.APpos]
%   This takes the same form as the spotsMat from flimscroll2 for
%   convenience. the 'nucleus' is arbitrary. No two spots are assigned to
%   the same nucleus in a single frame. jjust over time. cols with zeros
%   are meaningless
%
% Timothy Harden 2022

% initialize:
spotsCell = AllSpots.AllSpotsCells;
drift = DL.cumdriftlist;

% set up the output mat:
% number of frames:
frameN = DL.SequenceLentgh; % there is a typo in the script that makes DL
% frame range:
frameV = [DL.CorrectionRange(1):DL.CorrectionRange(2)];
% number of spots:
spotsN = 0;
for i = frameV
    n = spotsCell{i,2};
    spotsN = spotsN + n;
end
%allocate space:
spotsMat = zeros(spotsN,13);

% First, go frame by frame and populate the output matrix:
count = 1; % use this to keep track of spotsMat indicies
for i = frameV % loop through frames
    % get the number of spots in this frame:
    n = spotsCell{i,2};
    for j = 1:n % loop through spots in this frame
        % things that stay the same for each frame:
        % time (from the DL):
        spotsMat(count,1) = drift(i,4);
        % frame: 
        spotsMat(count,2) = spotsCell{i,3};
        % xpos:
        xpos = spotsCell{i,1}(j,1);
        xshift = drift(i,2);
        spotsMat(count,5) = xpos - xshift;
        % ypos:
        ypos = spotsCell{i,1}(j,2);
        yshift = drift(i,3);
        spotsMat(count,6) = ypos - yshift;
        % integreation stuff:
        spotsMat(count,7:10) = spotsCell{i,4}(j,:);
        % relative frame number (may need to change this): 
        spotsMat(count,11) = spotsCell{i,3};
        count = count + 1;
    end
end

thd = sqrt(thdx^2 + thdy^2); % this is the distance thingy
thd = 20; % for TSing
% assign numbers to spots, associate spots:
count = 1;
for i = 1:spotsN % loop through the entire length of spotsMat
% for i = 1:34
    if spotsMat(i,3) == 0
        % if there isn't a number assigned to this spot, assign one:
        spotsMat(i,3) = count;
        % list the distances between this spot and all subsequent:
        distList = [sqrt( (spotsMat(i + 1:end,5) - spotsMat(i,5)).^2 + (spotsMat(i + 1:end,6) - spotsMat(i,6)).^2 )];
        % are there any spots within the threshold distance?
        logi = distList < thd;
        % if so, assign them the same spot number:
        if sum(logi) > 0 
            % get the indicies of true values in logi:
            indV = find(logi);
            for j = 1:sum(logi)
                % get the indicie for this (subsequent) spot:
                ind = i + indV(j);
                spotsMat(ind,3) = count; % assign this spot the same number as spot i
            end
        end
        % update the spot assignment number:
        count = count + 1;
    end
end

% now put the spots back in their original spot for visual TSing in
% flimscroll:
count = 1;
for i = frameV % loop through frames
    % get the number of spots in this frame:
    n = spotsCell{i,2};
    for j = 1:n % loop through spots in this frame
        % things that stay the same for each frame:
        % xpos:
        xpos = spotsCell{i,1}(j,1);
        spotsMat(count,5) = xpos;
        % ypos:
        ypos = spotsCell{i,1}(j,2);
        spotsMat(count,6) = ypos;
        count = count + 1;
    end
end

mih = spotsMat;


