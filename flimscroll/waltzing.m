function out = waltzing(dat, frameNum, aoiNum, aoiinfo)
% out = waltzing(dat, frameNum, aoiNum, aoiinfo)
%
% This will shift AOIS in flimscroll for uploaded AOIS from tracking with
% gausian fitting. set FitChoice drop down to 'Moving AOI - from Tracking'
%
% INPUTS:
%   dat = aoifits.data from integrating aois with gauss fitting and
%      tracking. contained in handles.data in flimscroll
%   frameNum = ...
%   aoiNum = ...
%   aoiinfo = the aoiinfo list of AOI centers and initial frame number, as
%            contained in the handles.FitData of imscroll
%             [framenumber ave x y pixnum aoinumber]
%
% OUTPUTS:
%   out = a two member vecor, XYshift
%
% dat = [aoinumber framenumber amplitude xcenter ycenter sigma offset integrated_aoi (integrated pixnum) (original aoi#)]
%
% Harden 2018

%get the info for this aoi:
logi = dat(:,1) == aoiNum;
aoidat = dat(logi,:);

%get info for this frame
logi2 = aoidat(:,2) == frameNum;
framedat = aoidat(logi2,:);

%get the location of this aoi for this frame:
loc = framedat(4:5);

%original location:
logi3 = aoiinfo(:,6) == aoiNum;
oLoc = aoiinfo(logi3,3:4);

out = loc - oLoc;


