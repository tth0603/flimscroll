function mb = fitAspot(handles,aoiinfo)
% mb = fitAspot(handles,aoiinfo)
%
% A script to fit an individual spot, selected by the user, with a 2D gaussian
% in flimscroll under the addSpot callback. 
%
% INPUTS
%   handles: flimscroll gui handles
%   aoiinfo: an N x 6 matrix [frameNum numOfFramesAveraged xPos yPos pixelNum aoiNumber]
%
% OUTPUT
%   mb: a Nx6 matrix with cols: [1.gaussAmp 2.xCenter 3.yCenter 4.gaussSigma 5.offsetFromBackground 6.integratedIntensity];

% Timothy Harden 2019

% Fetch the current frame (appropriately averaged)
avefrm=getframes_v1Fli(handles);                       
% If the handles.BackgroundChoice is set to show the user
% a background-subtracted image, then use that background
% subtracted image in which to find spots.
if any(get(handles.BackgroundChoice,'Value')==[2 3])
     RollingBallRadius = str2double(get(handles.EditRollingBallRadius,'String'));
     RollingBallHeight = str2double(get(handles.EditRollingBallHeight,'String'));
     %Here to use rolling ball background (subtract off background)
     avefrm=avefrm-rolling_ball(avefrm,RollingBallRadius,RollingBallHeight);
elseif any(get(handles.BackgroundChoice,'Value')==[4 5])
     RollingBallRadius = str2double(get(handles.EditRollingBallRadius,'String'));
     RollingBallHeight = str2double(get(handles.EditRollingBallHeight,'String'));
     %Here to use Danny's newer background subtraction(subtract off background) 
     avefrm=avefrm-bkgd_image(avefrm,RollingBallRadius,RollingBallHeight);
end
 
spotNum = size(aoiinfo,1);
fits = [];
%loop through each added spot
for i = 1:spotNum
    %find the edges around the selected spot
    [xlowsml xhisml ylowsml yhisml]=AOI_Limits(aoiinfo(i,3:4),10); %this is a 20x20 AOI, should be plenty big but not too big
    %get just the area of the image around the spot center:
    firstaoi=avefrm(ylowsml:yhisml,xlowsml:xhisml); %weird naming convention from LJF
    %some params from ljf:
    mx=double( max(max(firstaoi)) );
    mn=double( mean(mean(firstaoi)) );
    %initial guess:
    inputarg0=[mx-mn 8 8 2 mn]; %%[amp centerx centery sigma offset]
    %gauss fit:
    gFit = gauss2dfit(double(firstaoi),inputarg0); %output: %[amp centerx centery sigma offset]; here the center is relative to the smaller image 'firstaoi'
    gFit = gFit';
    %now get the center of the spot in the ref frame of original image:
    %(top left corner of the small image [xlowsml ylowsml] plus the x & y
    %position of the spot center in the ref frame of the small image. 
    xCen = xlowsml + round(gFit(2)); %I dont like having to use the round feature here but not sure if its important and I am sure it doesn't really matters
    yCen = ylowsml + round(gFit(3));    
    %get the integrated intensity using the center of the spot as defined
    %by the fit (and integrating over the smaller 6x6 as you do in FindAllSpots.m:
    %get the spot limits:
    [xlowsml2 xhisml2 ylowsml2 yhisml2]=AOI_Limits(gFit(2:3),3);
    %harden 200116 to fit a sign error bug. kinda screws with spot position
    %but doesn't really effect anything
    if ylowsml2 < 0
        ylowsml2 = 0;
    end
    if xlowsml2 < 0
        xlowsml2 = 0;
    end
    %grab the pixel values:
    fitAOI=firstaoi(ylowsml2:yhisml2,xlowsml2:xhisml2); 
    %we use the smaller image 'firstaoi' to do the integration for convenience. 
    %This may have a pitfall if the point selected by the user is far from the 
    %actual spot center. To make up for this we draw a large box (20x20) around 
    %the user selected point in making 'firstaoi'.
    
    %sum the pixels to get integrated intensity
    integratedInt = sum(sum(fitAOI));
    %combine it all for putting out
    fits = [fits;gFit(1) xCen yCen gFit(4) gFit(5) integratedInt]; %[1.gaussAmp 2.xCenter 3.yCenter 4.gaussSigma 5.offsetFromBackground 6.integratedIntensity];
end
mb = fits;