function mh = AOIintegrate(handles)
% mh = AOIintegrate(handles)
%  
% OUTPUT
%   a mat of fit and integrated data for AOIs lacking a spot with cols: 
%   [1.time(s) 2.frameNum 3.nucleus 4.binarySpot 5.noSpotXpos 6.noSpotYpos 7.gaussAmp 8.sigma 9.offset 10.intInt 11.relativeFrameNumber]
%   that eventually gets save to analData in flimscroll2

% Timothy Harden 2020


% get the list of no spots:
noSpotMat = get(handles.loadNoSpots,'UserData'); %[1.time(s) 2.frameNum 3.nucleus 4.binarySpot 5.noSpotXpos 6.noSpotYpos 7.relativeFrameNumber]

% %get the image:
% im = get(handles.slider,'UserData');

FrameRange = min(noSpotMat(:,7)):max(noSpotMat(:,7)); %TH

[~, frmcol]=size(FrameRange);
%%%%%%%%%%%%%%%%%%%
currentimagenum=round(get(handles.slider,'Value')); %TH % Save the current image number value of the slider
noSpotFitMat = zeros(size(noSpotMat,1),11); %allocate space
noSpotFitMat(:,[1:6 11]) = noSpotMat(:,[1:7]); %make first six cols and last of output mat equal spot info

for frmindx=1:frmcol 
    if frmindx/100 == round(frmindx/10) %this is for reporting on progress
        frmindx
    end
    frame = FrameRange(frmindx);
    set(handles.slider,'Value',frame); % Set the slider value so that getframes_v1 will fetch the correct frame:
    guidata(gcbo,handles)  
    avefrm=getframes_v1Fli(handles); %grab the frame
    % get teh no spots in this frame
    frameLogi = noSpotMat(:,7) == frame;
    noSpotsInThisFrame = noSpotMat(frameLogi,:);
    aoirose = size(noSpotsInThisFrame,1);
    for ii = 1:aoirose
        %define the area around the aoi:
        [xlowsml xhisml ylowsml yhisml]=AOI_Limits(noSpotsInThisFrame(ii,5:6),5); %TH this is a 10x10 AOI, should be plenty big but not too big for gauss fit
        %get just the area of the image around the spot center:
        firstaoi=avefrm(ylowsml:yhisml,xlowsml:xhisml); %TH weird naming convention from LJF
        %some params from ljf:
        mx=double( max(max(firstaoi)) );
        mn=double( mean(mean(firstaoi)) );
        %initial guess:
        inputarg0=[mx-mn 5 5 2 mn]; %%[amp centerx centery sigma offset]
        %gauss fit:
        out = gauss2dfit(double(firstaoi),inputarg0); %TH output: %[amp centerx centery sigma offset]
        %to find the integreated intensity, use a much smaller aoi
        %size:
        [xlowsml2 xhisml2 ylowsml2 yhisml2]=AOI_Limits(noSpotsInThisFrame(ii,5:6),3);
        intAOI=avefrm(ylowsml2:yhisml2,xlowsml2:xhisml2);
        integratedInt = sum(sum(intAOI));
        %get the nuke number
        nuke = noSpotsInThisFrame(ii,3);
        %place the fit next to right nuke in the right frame in noSpotFitMat:
        noSpotFitMat(noSpotFitMat(:,11) == frame & noSpotFitMat(:,3) == nuke,7:10) = [out(1) out(4) out(5) integratedInt];
    end    
end
set(handles.slider,'Value',currentimagenum);                % Reinstate proper image number on slider

guidata(gcbo,handles)
mh=noSpotFitMat;   % [1.time(s) 2.frameNum 3.nucleus 4.binarySpot 5.noSpotXpos 6.noSpotYpos 7.gaussAmp 8.sigma 9.offset 10.intInt 11.relativeFrameNumber]
