function pc=FindAllSpotsFli(handles,maxnum)
%
%  function FindAllSpots(handles,maxnum)
% gfolder,frmrange,frmlimits,NoiseDiameter,SpotDiameter,SpotBrightness
% Will use spot picking algorithm to find all the spots in each frame
% defined by the frmrange.
% handles == handles structure of the calling program ( imscroll() )
% maxnum == maximum # of spots in each frame (for preallocating arrays:
%       too high means you have too many spots, nonspecific interactions or
%       some other problem)

% Copyright 2015 Larry Friedman, Brandeis University.
% edited by Harden 2018

% This is free software: you can redistribute it and/or modify it under the
% terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.

% This software is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
% A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with this software. If not, see <http://www.gnu.org/licenses/>.

%get inputs (TH)
nd = str2double(get(handles.noise,'String'));
b = str2double(get(handles.brightness,'String'));
%spot diam must be odd:
sd = round(str2double(get(handles.diameter,'String'))); 
if ~mod(sd,2)
    sd = sd+1;
end

% %get the image:
% im = get(handles.slider,'UserData');


%FrameRange=eval(get(handles.FrameRange,'String'));
FrameRange = eval([get(handles.FrameRange,'String') ';']); %TH
%FrameRange=ans;

[frmrose frmcol]=size(FrameRange);
                                            % preallocate the arrays for
                                            % storing spot information
% AllSpots.AllSpotsCells=cell(frmcol,3);  % Will be a cell array {N,3} N=# of frames computed, and
                                            % AllSpots{m,1}= [x y] list of spots, {m,2}= # of spots in list, {m,3}= frame#
AllSpots.AllSpotsCells=cell(frmcol,4);  %TH to include amp, sigma, offset from fitting as a mat in the 4th cell 
                                           
AllSpots.AllSpotsCellsDescription='{m,1}= [x y] list of spots in frm m, {m,2}= # of spots in list, {m,3}= frame#], {m,4} = [amp sigma offset integratedInt] list for spots in frame m';
AllSpots.FrameVector=FrameRange;         % Vector of frames whose spots are stored in AllSpotsCells
AllSpots.Parameters=[ nd sd b];  %TH
AllSpots.ParametersDescripton='[NoiseDiameter  SpotDiameter  SpotBrightness] used for picking spots';
AllSpots.aoiinfo2=get(handles.pick,'UserData');       % List of AOIs user has chosen      %TH
AllSpots.aoiinfo2Description='[frm#  ave  x  y  pixnum  aoi#]';



for indx=1:frmcol                   % Loop through all frame numbers
    AllSpots.AllSpotsCells{indx,1}(maxnum,:)=[0 0];   % Will contain the [x y ] list of spots, this preallocates space
    AllSpots.AllSpotsCells{indx,2}=0;                    % Will contain the number of spots for this frame
                                    % (nonzero entries of AllSpots{indx,1})
    AllSpots.AllSpotsCells{indx,3}=FrameRange(indx);   % Lists the frame number corresponding to these 3 cell arrays
end
% ave=1;  % Averaging number %TH
% pixnum=str2double(get(handles.PixelNumber,'String'));   % Pixel number

avefrm=getframes_v1Fli(handles);                       % Fetch the current frame(s) displayed
                                        % This is a dummy just to get size
[imagerose imagecol]=size(avefrm);                  % [ysize xsize]
xlow=1;xhigh=imagecol;ylow=1;yhigh=imagerose;         % Initialize frame limits
%TH: will need to do this if we magnify
if get(handles.showCrop,'Value')==1                  % Check whether the image magnified (restrct range for finding spots)  
    %limitsxy=eval( get(handles.MagRangeYX,'String') );  % Get the limits of the magnified region
    limitsxy = get(handles.aoiValue,'Value');                                               % [xlow xhi ylow yhi]
    xlow=limitsxy(1);xhigh=limitsxy(2);            % Define frame limits as those of 
    ylow=limitsxy(3);yhigh=limitsxy(4);            % the magnified region

end

currentimagenum=round(get(handles.slider,'Value')); %TH % Save the current image number value of the slider

for frmindx=1:frmcol            % Cycle through all frames, finding the spots
    if frmindx/500==round(frmindx/500) %this is for reporting on progress
        frmindx
    end
    %set(handles.get(handles.slider,'UserData');,'Value',FrameRange(frmindx));    % Set the slider value to current frame number so that getframes_v1 will fetch the correct frame
    set(handles.slider,'Value',FrameRange(frmindx));
 %TH
    guidata(gcbo,handles)
  
    avefrm=getframes_v1Fli(handles);                       % Fetch the current frame (appropriately averaged)
                % If the handles.BackgroundChoice is set to show the user
                % a background-subtracted image, then use that background
                % subtracted image in which to find spots.
%TH uncomment if you decide to mess with the background
     if any(get(handles.BackgroundChoice,'Value')==[2 3])
         RollingBallRadius = str2double(get(handles.EditRollingBallRadius,'String'));
         RollingBallHeight = str2double(get(handles.EditRollingBallHeight,'String'));
         % Here to use rolling ball background (subtract off background)            
         avefrm=avefrm-rolling_ball(avefrm,RollingBallRadius,RollingBallHeight);
     elseif any(get(handles.BackgroundChoice,'Value')==[4 5])
         RollingBallRadius = str2double(get(handles.EditRollingBallRadius,'String'));
         RollingBallHeight = str2double(get(handles.EditRollingBallHeight,'String'));
         % Here to use Danny's newer background subtraction(subtract off background) 
         avefrm=avefrm-bkgd_image(avefrm,RollingBallRadius,RollingBallHeight);
     end
    
    %dat=bpass(double(avefrm(ylow:yhigh,xlow:xhigh)),handles.NoiseDiameter,handles.SpotDiameter);
    dat=bpass(double(avefrm(ylow:yhigh,xlow:xhigh)),nd,sd);
    %pk=pkfnd(dat,handles.SpotBrightness,handles.SpotDiameter);
    pk=pkfnd(dat,b,sd);
    %pk=cntrd(dat,pk,handles.SpotDiameter+2);        % This is our list of spots in this frame FrameRange(frmindx)
    pk=cntrd(dat,pk,sd+2); %pk is [x y brightness (square of the radius of gyration))    
    
    [aoirose aoicol]=size(pk);
             
    if aoirose~=0       % If there are spots
        pk(:,1)=pk(:,1)+xlow-1;             % Correct coordinates for case where we used a magnified region
        pk(:,2)=pk(:,2)+ylow-1;
        %Th add gauss fitting here. make this optional in the GUI? 
        %allocate space:
        mb = zeros(aoirose,4);
        for ii = 1:aoirose
            %define the area around the aoi:
            [xlowsml xhisml ylowsml yhisml]=AOI_Limits(pk(ii,1:2),5); %TH this is a 10x10 AOI, should be plenty big but not too big
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
            [xlowsml2 xhisml2 ylowsml2 yhisml2]=AOI_Limits(pk(ii,1:2),3);
            fitAOI=avefrm(ylowsml2:yhisml2,xlowsml2:xhisml2);
            integratedInt = sum(sum(fitAOI));
            %asign to amp/sigma/offset/intInt to mat:
            mb(ii,:) = [out(1) out(4) out(5) integratedInt];
        end
        
        if aoirose>maxnum
            AllSpots.AllSpotsCells{frmindx,1}(1:maxnum,:)=pk(1:maxnum,1:2);   % List of spot centers
            AllSpots.AllSpotsCells{frmindx,2}=maxnum;                         % Number of detected spots we store
            AllSpots.AllSpotsCells{frmindx,3}=FrameRange(frmindx);            % Frame number 
            AllSpots.AllSpotsCells{frmindx,4}=mb; %[amp sigma offset]
        else
            AllSpots.AllSpotsCells{frmindx,1}(1:aoirose,:)=pk(1:aoirose,1:2);
            AllSpots.AllSpotsCells{frmindx,2}=aoirose;                         % Number of detected spots we store
            AllSpots.AllSpotsCells{frmindx,3}=FrameRange(frmindx);            % Frame number 
            AllSpots.AllSpotsCells{frmindx,4}=mb; %[amp sigma offset]
        end
    end
 
end

set(handles.slider,'Value',currentimagenum);                % Reinstate proper image number on slider

guidata(gcbo,handles)
pc=FreeAllSpotsMemoryFli(AllSpots);                        % Output structure containing cell array with spot record