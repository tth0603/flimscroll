function mih = rasterater4(varargin) 
% rasterater(analData1,...,analDataN,minLength,figNum) 
%
% This is for use with analData.nukesMat straight outta flimscroll. It will
% give a rastergram/heatmap output of single nucleus traces. Used to make a
% supplemental figure. 
%
% This differs from v3 in that it'll plot everything after gastrulation (or
% the end of a movie) as a different color
%
% INPUTS:
%   nukeMat == a nuke mat from the analData struc from flimscroll.
%       format: %[1.time(s) 2.frameNum 3.nucleus 4.binarySpot 5.nukeXpos 6.nukeYpos 7.spotGaussAmp 8.spotSigma 9.spotOffset 10.spotIntegratedIntensity 11.relativeFrameNumber 12.NC]
%   strtFrm == the start of the NC of interest. a scaler. THis is not
%   needed for data generated with v2 flimscroll and on
%   minLength == a scaler indicating the minimal length a trace must be to
%   be included (in frames)
%   figNum == where you gunna put it? (a scaler)
%   unsort == optional. any string. will not sort the nuclei intensity traces. ex:
%   'no sort'
%
% OUTPUTS
%   a raster plot of events (with heat map).
%
% USAGE
%   rasterater(analData.nukeMat,60,0,30)
%   rasterater(mat,60,0,30,'no sort')
%
% Harden 2018

% initiate:
repNum = length(varargin) - 2; % - 2 bc minFrame & figNum 
minLength = varargin{repNum + 1};
figNum = varargin{end};

% here is the color of the purgatory events, when reps end at different
% frames:
gastColr = 0.3;

% get teh relative final frame for each rep, the use the largest as the
% overall final frame:
finalFrameV = [];
nukeCell = {};
for j = 1:repNum 
    nukeCell{j} = varargin{j}.nukesMat;
    finalFrameV(j) = max(nukeCell{j}(:,2));
end
absFinalFrame = max(finalFrameV);

rasterMat = [];
for j = 1:repNum
    nukeMat = nukeCell{j};
    finalFrame = finalFrameV(j);
    % get the NaN's outta here 
    nukeMat(isnan(nukeMat)) = 0;
    % get those nukes:
    nukes = unique(nukeMat(:,3));
    for i = 1:length(nukes)
        mat2 = nukeMat(nukeMat(:,3) == nukes(i),:); %get all rose and cols for a given nuke
        if sum(mat2(:,4)) > 0 %any spots here? (of proper length)
            startPad = zeros(min(mat2(:,2)),1)'; %leading up to the first event 
            landingPad = zeros(finalFrame - max(mat2(:,2)),1)'; % between the last event and when we stop the analysis
            midPad = mat2(:,7)'; % the meat of it, the events
            gastPad  = ones(absFinalFrame - finalFrame,1)'.*gastColr;% an attempt at if the two reps end at a different time
            try 
                %rasterMat = [rasterMat;pad mat2(:,7)']; %cols become rose (plus we pad teh mat with zeros to account for all the frames in this NC where we arent makeing observations bc there are no spots
                rasterMat = [rasterMat;startPad midPad landingPad gastPad]; %cols become rose (plus we pad teh mat with zeros to account for all the frames in this NC where we arent makeing observations bc there are no spots
            catch
                keyboard
            end
        end
    end
end

%attempt to limit frame length (spoiler: it works.)
rose = size(rasterMat,1);
if minLength > 0
    mb = [];
    for j = 1:rose
        %here is the events vector for nucleus v
        v = rasterMat(j,:);
        %find the runs of events, where they start and stop
        notZero = [false, (v > 0),false];
        indx = [strfind(notZero, [false, true]); ...
                strfind(notZero, [true, false]) - 1]; %a 2 x noe sized mat
        %num of events:
        noe = size(indx,2);
        %length of events (a vector of length noe)
        el = indx(2,:) - indx(1,:) + 1;
        for m = 1:noe
            %if an event for nucleus j is shorter than minLength, make it
            %zero:
            if el(m) < minLength
                v(indx(1,m):indx(2,m)) = 0;
            end
        end
        %only put the events vector back into rasterMat if there are still
        %events described in it:
        if sum(v) > 0 
            mb = [mb;v];
        end
    end
    rasterMat = mb;
end                
    
%sort rasterMat (there are a million ways to do this. this shitty soln prolly isnt the best):
indV = [];
for k = 1:size(rasterMat,1)
    [~,col] = (find(rasterMat(k,:),1)); %find first col with a non zero value in each row
    indV = [indV;col]; %place these column indices in a vector the length of rasterMat
end
%append these col indices to the first column of rasterMat
rasterMatInd = [indV rasterMat];
% sort this based on teh first col:
sortMat = sortrows(rasterMatInd,1);
%get rid of the first col (the indicies)
rasterMatSort = sortMat(:,2:end);

% if you dont want to have a sorted output, simply edit to reassign rasterMatSort
% to rasterMat

%plot
%get a time vector to replace the frames on x axis
%ts = unique(nukeMat(:,1)) - nukeMat(1,1); %subtract off initial time, as we go from zero at start of nc
tStep = mean(diff(unique(nukeMat(:,1))))*2; %so many parens. this is the average time btwn acqns
% mxT = (max(nukeMat(:,11)) + size(pad,2))*tStep; %the is the time of the last frame is s relative to the start of the nc
mxT = finalFrame*tStep;
ticStep = 20; %every tenth time point place a tick mark on the x axis, can change this
%plotTs = round(ts(1:ticStep:end)); 
plotTs = round(0:tStep*ticStep:mxT);
% tics = 0:ticStep:max(nukeMat(:,11)) + size(pad,2); %for the actual tick placement, go every 10 frames starting at zero and ending at the max number of (relative) frames PLUS the number of frames we've padded the traces with to account for the discrepency of nc start frame and the first frame tracked
tics = 0:ticStep:finalFrame;

%%%%%% make the mask binary:
% use linear indexing here:
[x,y] = size(rasterMatSort);
for i = 1:x*y
    if rasterMatSort(i) > 0 && rasterMatSort(i) ~= gastColr
        rasterMatSort(i) = 1;
    end
end

%%%%%% make the area prior to first passage gray:
kolor = 0.8; % the value that we want to make the first passage area 
% rose = size(rasterMat,1);
for i = 1:rose %cycle through rose bc I don't know how to do this for the entire matrix
    ind = find(rasterMatSort(i,:),1,'first');  %linear indexing to find the first event in row j
    ind = ind - 1; %we want all the indicies prior to the first event
    rasterMatSort(i,1:ind) = kolor;
end

%%%% make the area after gastrulation dark grey:



figure(figNum);
imagesc(rasterMatSort);shg %axis('equal');shg
j = gray; 
colormap(j)
set(gca,'XTick',tics,'xticklabel',plotTs); %,'XTickLabelRotation',45);
end

