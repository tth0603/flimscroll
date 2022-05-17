function rasterater(nukeMat,strtFrm,minLength,figNum,varargin) 
% rasterater(nukeMat,strtFrm,minLength,figNum,'no sort') 
%
% This is for use with analData.nukesMat straight outta flimscroll. It will
% give a rastergram/heatmap output of single nucleus traces.
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

% %if the start frame is already contained in the nukeMat
% %make a zero vector the length of the number of frames btwn nc14 and the
% %first frame of observation:
% if isempty(strtFrm)
%     pad = zeros(min(nukeMat(:,2)),1)';
% else
%     pad = zeros(min(nukeMat(:,2)) - strtFrm,1)';
% end

% get the max length of observation
finalFrame = max(nukeMat(:,2));
% finalFrame = 149; %edit this to make x axis a certain size

% get the NaN's outta here 
nukeMat(isnan(nukeMat)) = 0;
% get those nukes:
nukes = unique(nukeMat(:,3));
rasterMat = [];
for i = 1:length(nukes)
    mat2 = nukeMat(nukeMat(:,3) == nukes(i),:); %get all rose and cols for a given nuke
    if sum(mat2(:,4)) > 0 %any spots here? (of proper length)
        startPad = zeros(min(mat2(:,2)),1)';
        landingPad = zeros(finalFrame - max(mat2(:,2)),1)';
        try 
            %rasterMat = [rasterMat;pad mat2(:,7)']; %cols become rose (plus we pad teh mat with zeros to account for all the frames in this NC where we arent makeing observations bc there are no spots
            rasterMat = [rasterMat;startPad mat2(:,7)' landingPad]; %cols become rose (plus we pad teh mat with zeros to account for all the frames in this NC where we arent makeing observations bc there are no spots
        catch
            keyboard
        end
    end
end

%attempt to limit frame length (spoiler: it works.)
if minLength > 0
    mb = [];
    for j = 1:size(rasterMat,1)
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
%if you want unsorted output:
if exist('varargin{1}','var') && ischar(varargin{1})
    rasterMatSort = rasterMat;
end
%plot
%get a time vector to replace the frames on x axis
%ts = unique(nukeMat(:,1)) - nukeMat(1,1); %subtract off initial time, as we go from zero at start of nc
tStep = mean(diff(unique(nukeMat(:,1)))); %so many parens. this is the average time btwn acqns
% mxT = (max(nukeMat(:,11)) + size(pad,2))*tStep; %the is the time of the last frame is s relative to the start of the nc
mxT = finalFrame*tStep;
ticStep = 10; %every tenth time point place a tick mark on the x axis, can change this
%plotTs = round(ts(1:ticStep:end)); 
plotTs = round(0:tStep*ticStep:mxT);
% tics = 0:ticStep:max(nukeMat(:,11)) + size(pad,2); %for the actual tick placement, go every 10 frames starting at zero and ending at the max number of (relative) frames PLUS the number of frames we've padded the traces with to account for the discrepency of nc start frame and the first frame tracked
tics = 0:ticStep:finalFrame;

figure(figNum);
imagesc(rasterMatSort,[0 300]);shg %axis('equal');shg
j = parula; 
%make bkgrnd white
j(1,:) = [1 1 1];
colormap(j)
colorbar;
ylabel('nucleus number');
xlabel('time (s)');
set(gca,'XTick',tics,'xticklabel',plotTs,'XTickLabelRotation',45);

end

