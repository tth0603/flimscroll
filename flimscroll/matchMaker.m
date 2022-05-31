function mib = matchMaker(aoifits,AllSpots,varargin)
% out == matchMaker(aoifits,AllSpots,timeBase,ncs)
%
% This script is called in flimscroll at the "Pair spots" button press.
% Takes in the aoifits and AllSpots strucs from nuclei fitting (Tracks botton)
% and spot picking (Find all spots button), respectively. 
% 
% The output defines the struc upon which all DS analysis relies. This is
% called in flimscroll by handles.makeMatches,'UserData'. and saved in imscroll
% as 
% this is still a 
% The
% content and format of this is still a work in progress. 
%
%want: [frame time nucleus binarySpot spotAmplitude spotSigma spotOffset spotIntensity xLocalLocation yLocalLocation percentAP] 
% presently getting: [time frame nucleus binarySpot distanceFromNucleus xPos yPos amplitude sigma offset intensity relativeFrameNumber]
%Harden 2018

%get info:
% aoifits.dataDescription
% [aoinumber framenumber amplitude xcenter ycenter sigma offset integrated_aoi (integrated pixnum) (original aoi#)]
% AllSpots.AllSpotsCellDescription
% '{m,1}= [x y] list of spots in frm m, {m,2}= # of spots in list, {m,3}= frame#], {m,4} = [amp sigma offset] list for spots in frame m';
spotsMatDescription = '[1.timeSinceDivision(s) 2.frameNumSinceDivision 3.nucleus 4.distanceFromNucleusCntr 5.xPos 6.yPos 7.gaussAmp 8.sigma 9.offset 10.integratedIntensity 11.absFrameNumber 12.NC]';
nukesMatDescription = '[1.timeSinceDivision(s) 2.frameNumSinceDivision 3.nucleus 4.binarySpot 5.nukeXpos 6.nukeYpos 7.spotGaussAmp 8.spotSigma 9.spotOffset 10.spotIntegratedIntensity 11.absFrameNumber 12.NC]';

spots = AllSpots.AllSpotsCells;
nukes = aoifits.data;
frameNum  = length(spots); % number of frames

%get the frame range over which this stuff has been integrated (many ways
%to do this):
startFrm = spots{1,3};
endFrm = spots{end,3};

%if there is no time base, put in a vector of zeros (not sure this is
%relevant any longer, so I commented it out)
% if length(varargin) == 1
    ttb = varargin{1};
% else 
%     ttb = zeros(frameNum,1);
% end
timeBase = ttb(startFrm:endFrm);

%grab NC input
ncs = varargin{2};
%get rid of nans:
ncs = ncs(~isnan(ncs(:,2)),:);
%how many cycles do we have data for?
ncNum = size(ncs,1);
%tack another row so we can slip the last NC in:
ncs = [ncs;1 length(ttb)+1]; %ignore the 1. thats jsut a place holder. the '+1' avoids fence post error
%make a mat 2 x frame num long with the nc for each frame:
ncMat = [];
%first take care of the frames running up to the first NC with an input?
if ncs(1,2) > 1
    frms1 = [1:ncs(1,2)-1]';
    frmLen1 = length(frms1);
    currNc1 = ones(frmLen1,1)*(ncs(1,1) - 1);
    ncMat = [currNc1 frms1];
end
%then the rest of the NCs
for i = 1:ncNum
    frms = [ncs(i,2):ncs(i+1,2)-1]'; %this indeving only works bc we slipped in the extra row above
    frmLen = length(frms);
    currNc = ones(frmLen,1)*(ncs(i,1));
    ncMat = [ncMat;currNc frms];
end
%since we cannot presently analyze across multiple ncs, lets figure
%relative times and frames here:
nc = ncMat(ncMat(:,2) == spots{1,3},1); %here we use logi indexing to pluck out the NC from ncMat using the fact that spots{i,3} is the absolute frame number
%get the starting frame of the corresponding NC:
zeroFrame = ncs(ncs(:,1) == nc,2);
%bc nc & ttb have the same length, we get the abs start time of the NC via:
zeroTime = ttb(zeroFrame);

%heres a mat where every spot gets a row. we likely also want one where
%every nucleus gets a row, for fraction of nuclei, etc...
spotsMat = []; %so far, this runs fast so no need to allocate space
altMat = [];
for i = 1:frameNum  %indexing for each frame
    frameSpot = spots{i,1}; %get spot cell for this frame
    nc = ncMat(ncMat(:,2) == spots{i,3},1); %here we use logi indexing to pluck out the NC from ncMat using the fact that spots{i,3} is the absolute frame number
    if ~isempty(frameSpot) %look for spots in this frame
        frameSpotNum = spots{i,2}; %num of spots in this frame
        frameLogi = nukes(:,2) == spots{i,3}; %3rd col of spots cell array is frame number
        frameNukes = nukes(frameLogi,:); %get all nuclei info in this frame, this should be the total num of nukes for this motion picture
        spotFitParams = spots{i,4}; %this is the new bit to AllSpots cell array that flimscroll adds in FindAllSpotsFli
        for j = 1:frameSpotNum %indexing for each spot in the frame
            spotLoc = frameSpot(j,:); %each spot
            spotFit = spotFitParams(j,:); %[amp sigma offset integratedInt]
            disList = [sqrt((frameNukes(:,4) - spotLoc(1)).^2 + (frameNukes(:,5) - spotLoc(2)).^2) frameNukes(:,1)]; %[distance nukNumber] %only use sqrt for TSing
            %above computes the distance between the spot center and each nuclei in the frame
            [sortdistance I]=sort(disList(:,1));   % Now sort list (ascending)
            nukMatch = disList(I(1),:); %here we get the distance and nuk number ; outside TSing, this can just be disList(I(1),:);
            spotsMat = [spotsMat;[timeBase(i)-zeroTime spots{i,3}-zeroFrame nukMatch(2) nukMatch(1) spotLoc(1) spotLoc(2) spotFit(1) spotFit(2) spotFit(3) spotFit(4) spots{i,3} nc]]; %added the distance (nukMatch(1)) for TSing%[1.time 2.frameNum 3.nucleus 4.distanceFromNucleusCntr 5.xPos 6.yPos 7.gaussAmp 8.sigma 9.offset 10.integratedIntensity 11.relativeFrameNumber]
            altMat = [altMat;[timeBase(i)-zeroTime spots{i,3}-zeroFrame disList(I(2),2) disList(I(2),1) disList(I(3),2) disList(I(3),1) disList(I(4),2) disList(I(4),1)]]; %This is to get the 2nd closest nuke (and 3rd), for cases when multiple spots are assigned to the same nuke
        end
    end
end
%let's take care of multiple spots assigned to the same nucleus by keeping
%the nearest spot and assigning the other(s) to the next nearest neighbor
%for the case of 2 spots/nuke:
for i = 1:length(spotsMat)
    nuk = spotsMat(i,3);
    frm = spotsMat(i,2);
    misMat = spotsMat(spotsMat(:,3) == nuk & spotsMat(:,2) == frm,:);
    if size(misMat,1) > 1 %if we find more than one spot at that nuke at that time
        nearest = min(misMat(:,4)); %find the closest distance
        if nearest ~= spotsMat(i,4)
            spotsMat(i,[3 4]) = altMat(i,[3 4]); % if the current spot isnt the closest, assign that spot to the next nearest neighbor
        end
    end
end
%i've tried to iterate this to fix it when this attempt fails. Thusfar that
%has not worked and I am not ATM clever enough to come up with a better
%optimization. My fix now is to allow the user to correct these via the
%Correct spot assignment panel
    
%here's a mat including all nuclei.
%I added this later and cant read code, so that's why they are written funny
nukeNum = size(aoifits.centers,1);
nukeLen = size(nukes,1);
%nukes mat is the same length as aoifits.data, with 12 cols. we go ahead
%and allocate that
nukesMat = zeros(nukeLen,12);
% time:
% nukesMat will be first indexed by nucleus, then time/frame. make a
% numOfNukes x numOfTimePoints sized mat. then linearize it and shove it
% into nukesMat in a subsequent for loop
nukeTime = [];
for i = 1:max(nukes(:,1)) %1:number of nuclei in track
    nukeTime = [nukeTime;ttb(startFrm:endFrm)']; %a numOfNukes x numOfTimePoints mat
end
nukeTime = nukeTime - zeroTime; %a sloppy correction due to post hoc edits makes the times relative to the start of the NC. 
%same for the relative frame:
absFrame = [];
for i = 1:max(nukes(:,1))
    absFrame = [absFrame;startFrm:endFrm];
end
%rip some stuff from aoifits.data, rearrange to preseve consistent notation
%with spotsMat:
nukesMat(:,[2 3 5 6]) = nukes(:,[2 1 4 5]);
%Make the frame realtive to the start of the relevan NC:
nukesMat(:,2) = nukesMat(:,2) - zeroFrame;
%now the binary bit, along with spot characteristics if we have them:
for i = 1:nukeLen
    %assign time:
    nukesMat(i,1) = nukeTime(i); %linear indexing on nukeTime works bc of the way we made the mat above
    %assign relative frame:
    nukesMat(i,11) = absFrame(i); %linear indexing
    logi1 = nukesMat(i,1) == spotsMat(:,1) & nukesMat(i,3) == spotsMat(:,3); %get the matching frame (col 2) and nucleus (col 3)
    bin = sum(logi1); %sum the logi vector. this should give a scaler of either 1 or 0. if >1, see below
    %include binary
    nukesMat(i,4) = bin;
    %and spot characteristics:
    if bin == 1 %if we find a spot for that nuke at that time
        mb = spotsMat(logi1,:);
        nukesMat(i,[7 8 9 10]) = mb(:,[7 8 9 10]);
    elseif bin == 0 %if we do not find a spot
        nukesMat(i,[7 8 9 10]) = NaN;
    elseif bin > 1 %if multiple spots get assigned to a nuke in a single frame we must figure that out what to do
        %update, this should be taken care of 
        mb = spotsMat(logi1,:);
        %currently we assign the closest spot:
        logi2 = sum(mb(:,4:5).^2,2) == min(sum(mb(:,4:5).^2,2)); %this'll pick out the closest spot
        try
            nukesMat(i,[7 8 9 10]) = mb(logi2,[7 8 9 10]);
        catch
            keyboard
        end
        %nukesMat(i,[7 8 9 10]) = mb(mb(:,4) == min(mb(:,4)),[7 8 9 10]);   %logical indexing without bothering to make a logical object
        %warn the user:
        nuke = nukesMat(i,3);
        frame = nukesMat(i,2);
        fprintf('multiple spots were found for nuclues %d at frame %d \n',nuke,frame);
    end
    %fianlly, add nc:
    nc = ncMat(ncMat(:,2) == nukesMat(i,2),1);
    nukesMat(i,12) = nc;
end
fprintf('for instances when multiple spots were found for a single nucleus, \n we assign the closest spot to that nucleus in analData.nukesMat. \n However there will be spots in analData.spotsMat that are assigned to \n the same nuke in the same frame. These can be dealt with post processing. \n')

%add intervals
intsStruc = makeInts(nukesMat);

%assign outputs
mib.spotsMat = spotsMat;
mib.spotsMatDescription = spotsMatDescription;
mib.nukesMat = nukesMat;
mib.nukesMatDescription = nukesMatDescription;
mib.cia = intsStruc.cia;
mib.cumulativeIntervalsArrayDescription = intsStruc.cumulativeIntervalsArrayDescription;



%some notes: 
% AllSpots has a cell array that presently contains some of what we want, less
% any fitting (that is each spot and its loction in each frame). this is 
% the cell array AllSPots.AllSpotsCell.
% the cell array format may be unfortunate, but its probably
% the best format. we then must match that to the nuclei postitions from
% aoifits.data structure. this is an easier format but is probably not the
% most efficient way.