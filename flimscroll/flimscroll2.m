function varargout = flimscroll2(varargin)
% FLIMSCROLL2 MATLAB code for flimscroll2.fig
%      FLIMSCROLL2, by itself, creates a new FLIMSCROLL2 or raises the existing
%      singleton*.
%
%      H = FLIMSCROLL2 returns the handle to a new FLIMSCROLL2 or the handle to
%      the existing singleton*.
%
%      FLIMSCROLL2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FLIMSCROLL2.M with the given input arguments.
%
%      FLIMSCROLL2('Property','Value',...) creates a new FLIMSCROLL2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before flimscroll2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to flimscroll2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% 
% Where data is stored within the GUI:
%       spots/nuclei from feature picker: handles.pick,'UserData' (aoiinfo2 mat)
%       time base handles.timebase,'UserData' (a struc named timebaseStruc containing a timebase vector ttb)
%       nuclei locations: handles.pick,'UserData' (aoiinfo2 mat)  
%       nuclei tracking: handles.trackNuclei,'UserData' (aoifits struc)
%       spots: handles.findAllSpots,'UserData' (AllSpots struc)
%       spot-nuclei pairs: handles.makeMatches,'UserData' (this is presently a struc called analData with a spotsMat & a nukesMat along with their descriptions. made in script 'matchMaker' and added to by the 'Add AP pos.' button)
%       file locations: handles.outputText,'UserData' (set at GUI launch)
%       maxProj directory: handles.ms2dir,'UserData',{fp fn}
%       BF readers: handles.load,'UserData' (MS2) handles.load2,'UserData' (his)
%       getting input frame range: frms = eval([get(handles.FrameRange,'String') ';']); (a vector of length number of frames)
%       finding the total number of frames: get(handles.slider,'Max') (a scaler)
%       driftlist: handles.createDriftlist,'UserData',DL (this is to improve nuclei treacking.)
%       AP position: handles.findAP,'UserData' (a struc with two members: coordA = [xpos ypos] & coordP = [xpos ypos])
%       analData: handles.makeMatches,'UserData' (a struc that contains all the data to ananlyze. presently, it contains the spotsMat & nukesMat matrices that contain all quantified information for the MS2 spots and nuclei, respectively. it also contains Descreiption strings for each)
%       Integrated intensity: This is calculated in two locations:
%           findAllSpots (via FindAllSpotsFli.m) and addSpot (via fitAspot.m).
%           Presently we integrate over a 6x6 AOI centered around the gaussian
%           peak of the spot.
%
%   Bodies buried:
%       1. cannot track through mitosis. each n.c. must be analyzed
%       individually. --will likely never be addressed
%       3. Only can handle max proj data. will need to be improved. 
%       7. get rid of nuclei doublets
%       
% Timothy Harden 2018
% 

% Last Modified by GUIDE v2.5 27-Aug-2021 14:01:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @flimscroll2_OpeningFcn, ...
                   'gui_OutputFcn',  @flimscroll2_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before flimscroll2 is made visible.
function flimscroll2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no outputText args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to flimscroll2 (see VARARGIN)

% Choose default command line outputText for flimscroll2
handles.outputText = hObject;

%just initially, this will change when images are loaded:
set(handles.slider,'Max',1000);
set(handles.slider,'Min',1);
set(handles.slider,'Value',1);
set(handles.slider,'SliderStep',[1/(1000-1) 10/(1000-1)])

%set file locations
load('fliFileLocations.dat','-mat','fileLocations')
set(handles.outputText,'UserData',fileLocations);

%turn off axes:
axes(handles.axes1);
axis('off');

%initialize some stuff that is not is handles. esp good for fitting
% handles.FitData = [];
% handles.EditRollingBallRadius = 15;
% handles.EditRollingBallHeight = 5; %if Choice is a thing, then we can also switch these
handles.Pixnums = [];
% handles.TrackAOIs = 1; %this is now the tag for the hidden text titled "Fit Choice"
handles.TiffFolder = '';

version = 2;
fprintf('version %d \n',version);
%v2: 1. added "Correct spot assigment" panel
% 2. within analData time and frame is now relative to start of NC; abs
% frame # is now col 11.
% 3. added correct-a-track panel
% 4. magnify view can now be set by manual input
% 5. can now add/remove spots 
% 6. background subtraction enabled for disdplay and spot picking
% 7. Added two condition spot picking


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes flimscroll2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = flimscroll2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning outputText args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line outputText from handles structure
varargout{1} = handles.outputText;

% --- Executes on slider movement.
function slider_Callback(hObject, eventdata, handles)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
%grab slider value
imagenum = round(get(hObject,'Value'));
display = num2str(imagenum);
%set displayed frame number
set(handles.frameNum,'String',display)
set(handles.frameNum,'Value',imagenum)

%crop stuff
crp = get(handles.showCrop,'Value');
% aoi = get(handles.aoiValue,'Value'); %old. try this if mag breaks
if crp == 1
    aoi = eval([get(handles.aoiValue,'String') ';']); %now [xmin xmax ymin ymax]
end

%decide if it's MS2 or his:
hisVal = get(handles.hisImage,'Value');
ms2Val = get(handles.ms2image,'Value');

%show numbers?
shoNum = get(handles.aoiNum,'Value');

%feature display size
pixnum = str2double(get(handles.PixelNumber,'String')); %how big of spots

%get the channel:
if ms2Val
    %get bf reader:
    reader = get(handles.load,'UserData');
elseif hisVal
    reader = get(handles.load2,'UserData');
else 
    set(handles.instrucs,'String','something went wrong with MS2/his toggle')
end

%get the image:
%im = bfGetPlane(reader,reader.getIndex(0, 1, imagenum - 1));  %in the past this has worked. I am growing a bit tired of how finicky this bf shit is
im = bfGetPlane(reader,reader.getIndex(imagenum - 1, 0, 0)+1);       %getIndex is (Z, C, T) I think. the channel should be 0 for gfp, but may change
%and check if you'd like to remove background stuff
%get background subtraction stuff:
if get(handles.BackgroundChoice,'Value') ~= 1
    RollingBallRadius = str2double(get(handles.EditRollingBallRadius,'String'));
    RollingBallHeight = str2double(get(handles.EditRollingBallHeight,'String'));
end
switch get(handles.BackgroundChoice,'Value')
    case 1 %raw imagw
    case 2 %im og the background to be subreacted
        im=rolling_ball(im,RollingBallRadius,RollingBallHeight);
    case 3 %im after background subtraction
        im=double(im)-double(rolling_ball(im,RollingBallRadius,RollingBallHeight));
    case 4 %method 2 
        im=bkgd_image(im,RollingBallRadius,RollingBallHeight);
    case 5
        tic;
        im=double(im)-double(bkgd_image(im,RollingBallRadius,RollingBallHeight));
        toc
end

%get the image scale stuff
if get(handles.autoScale,'Value') == 0
    clowval = get(handles.minSlider,'Value');
    chival = get(handles.maxSlider,'Value');
else
    clowval = min(min(im));
    chival = max(max(im));
    %make the display informative:
    lowStr = round(clowval);
    hiStr = round(chival);
    %so the max scale slider behaves nicely with high chival:
    mxx = get(handles.maxSlider,'Max');
    if chival <= mxx
        set(handles.maxSlider,'Value',chival);
        set(handles.maxScale,'String',hiStr);
    else
        set(handles.maxSlider,'Value',mxx);
        set(handles.maxScale,'String',hiStr);
    end
    if clowval < 0 
        set(handles.minSlider,'Value',0)
        set(handles.minScale,'String',0);
    else 
        set(handles.minSlider,'Value',clowval); 
        set(handles.minScale,'String',lowStr);
    end      
end
%playing around:
% im = imadjust(im,[0; 1],[0; 1],2);

%plot/display
hAx = handles.axes1;
if crp && length(aoi) == 4
    %smlIm = im(aoi(2):aoi(4)+aoi(2),aoi(1):aoi(3)+aoi(1)); %weird format of aoi bc the axis of the image is inverted, i think
    %imshow(imadjust(mat2gray(smlIm,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);
    %imagesc(im,[clowval chival] );axis('equal');colormap(gray(256));axis([aoi(1),aoi(1)+aoi(3),aoi(2),aoi(2)+aoi(4)])  %axis([xmin xmax ymin ymax])
    %old cropping:
    %imagesc(im,[clowval chival] );axis('equal');axis('off');colormap(gray(256));axis([aoi(1),aoi(1)+aoi(3),aoi(2),aoi(2)+aoi(4)])  %axis([xmin xmax ymin ymax])
    %new cropping 
    imagesc(im,[clowval chival] );axis('equal');axis('off');colormap(gray(256));axis(aoi)
else
%         imshow(imadjust(mat2gray(im,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);
    imagesc(im,[clowval chival] );axis('equal');axis('off');%colormap(gray(256));%axis(limitsxy)
    %imagesc(im,[0 256] );axis('equal'); colormap(gray(256));%axis(limitsxy)

end
set(handles.slider,'UserData',im);

%show tracked nuclei:
if get(handles.displayTracks,'Value') == 1
    aoifits = get(handles.trackNuclei,'UserData');
    %analData = get(handles.makeMatches,'UserData'); ...%incomplete. an
    %attempt to avoid loading tracks snd just get them from analData
    if ~isempty(aoifits)
        aoirose = size(aoifits.aoiinfo2,1);
        data = aoifits.data;
        aoiinfo2 = aoifits.aoiinfo2;
        %make sure we have data for the current frame:
        fitMx = max(data(:,2));fitMn = min(data(:,2));
        if imagenum >= fitMn && imagenum <=fitMx
            for indx=1:aoirose
                XYshift = waltzing(data, imagenum, indx, aoiinfo2); %inputs: a aoifits.data mat; the frame number; the aoi number; aoiinfo2 mat
                viscircles(aoiinfo2(indx,3:4)+XYshift,pixnum/2,'Color','y','LineWidth',0.5,'EnhanceVisibility',false);
                if shoNum == 1 %do you want to put the AOI number on the image?
                    text(aoiinfo2(indx,3) + XYshift(1),aoiinfo2(indx,4) + XYshift(2),num2str(aoiinfo2(indx,6)),'FontSize',8,'Color','y');
                end
            end  
        end
    %elseif ~isempty(analData)
        %nukesMat = analData.nukesMat;
        %... incomplete
    end
end
%show nuclei/or features from spot picker:
if get(handles.displayNuclei,'Value') == 1
    FitData = get(handles.pick,'UserData');
    %make sure we ahve features to show
    if ~isempty(FitData)
        aoirose = size(FitData,1);
        if get(handles.displayDriftlist,'Value') == 0
            for indx=1:aoirose %indexing on aoi number
                viscircles(FitData(indx,3:4),pixnum/2,'Color','b','LineWidth',0.5,'EnhanceVisibility',false);
                if shoNum == 1
                    text(FitData(indx,3),FitData(indx,4),num2str(FitData(indx,6)),'FontSize',8,'Color','y');
                end
            end 
        else %adjust positions if you want to track drift:
            for indx=1:aoirose
                DL = get(handles.createDriftlist,'UserData');
%                 shift = ShiftAOIFli(FitData(1,6),indx,FitData,DL);
                shift = ShiftAOIFli(indx,imagenum,FitData,DL); %(aoiNum frameNum aoiinfo2 driftlist)
%                 aoiinf(logic,3:4)=aoiinf(logic,3:4)+ShiftAOIFli(aoiinf(1,6),indxx,FitData,DL);
                viscircles(FitData(indx,3:4)+shift,pixnum/2,'Color','b','LineWidth',0.5,'EnhanceVisibility',false);
                if shoNum == 1
                    text(FitData(indx,3)+shift,FitData(indx,4),num2str(FitData(indx,6)),'FontSize',8,'Color','y');
                end
            end 
        end
    end
end

%see if there are spots (matched or otherwise) and if we want to plot them:
if get(handles.plotSpots,'Value') == 1
    %for paired spots
    analData = get(handles.makeMatches,'UserData');  %this is where the matched spots struc lives. current form: [frame nucleus binarySpot distance xPos yPos]
    %for regular spots
    spots = get(handles.findAllSpots,'UserData');
    
    if ~isempty(analData)
        %this is where the matched spots struc lives. current form: [1.relativeTime 2.relativeFrameNum 3.nucleus 4.distanceFromNucleus 5.xPos 6.yPos 7.amplitude 8.sigma 9.offset 10.intensity 11.relativeFrameNumber 12.APpos]
        spotPairs = analData.spotsMat;
        spotsLogi = spotPairs(:,11) == imagenum;
        currentSpots = spotPairs(spotsLogi,:); %get spots corresponding to this frame
        frameSpotsLen = size(currentSpots,1);  %get number of these spots
        axes(handles.axes1);
        hold on;
        for i = 1:frameSpotsLen
            viscircles(currentSpots(i,5:6),3,'Color','r','LineWidth',0.5,'EnhanceVisibility',false);
            if shoNum == 1
                text(currentSpots(i,5)+5,currentSpots(i,6)+3,num2str(currentSpots(i,3)),'FontSize',8,'Color','r'); %x & Y shifted by 3 pix, write out the aoi number
            end
        end
        hold off
    elseif ~isempty(spots)  %place plotting regular spots here
        spotsCell = spots.AllSpotsCells;
        spotsCellLen = size(spotsCell,1);
        for i = 1:spotsCellLen
            if imagenum == spotsCell{i,3}
                radVec = ones(size(spotsCell{i,1},1),1).*3; %make them radius 3
                viscircles(spotsCell{i,1},radVec,'Color','r','LineWidth',0.5,'EnhanceVisibility',false); %take advantage of the fact that you can input mats into viscircles
            end
        end
    end
end

% if there are no spots
if get(handles.displayNoSpots,'Value') == 1
   %check analData first (this isn't set up yet)
   
   %otherwise use noSpots loaded in Spots panel
   noSpotMat = get(handles.loadNoSpots,'UserData'); % [1.time(s) 2.frameNum 3.nucleus 4.binarySpot 5.noSpotXpos 6.noSpotYpos 7.relativeFrameNumber]
   % get the no spots for this frame:
   noSpotsLogi = noSpotMat(:,7) == imagenum;
   currentNoSpots = noSpotMat(noSpotsLogi,:);
   currentNoSpotsLen = size(currentNoSpots,1);
   axes(handles.axes1);
   hold on;
   for i = 1:currentNoSpotsLen
       viscircles(currentNoSpots(i,5:6),3,'Color',[0.929 0.694 0.125],'LineWidth',0.5,'EnhanceVisibility',false); %make them orange
       if shoNum == 1
            text(currentNoSpots(i,5)+5,currentNoSpots(i,6)+3,num2str(currentNoSpots(i,3)),'FontSize',8,'Color',[0.929 0.694 0.125]); %x & Y shifted by 3 pix, write out the aoi number
       end
   end
end

% --- Executes during object creation, after setting all properties.
function slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)  %loading for MS2
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%in case the slider has been fucked with:
set(handles.slider,'Value',1);

currDir = cd;
%if a file has already been retrieved, go to that dir
try
    dr = get(handles.hisDir,'UserData');
    cd(dr{1});
end

%get the image file
[fn fp] = uigetfile('*.tiff','pick your MS2 max proj tif');
reader = bfGetReader([fp fn]);
im = bfGetPlane(reader, 1); %1st image

%plot 1st image to confirm load
ms2Val = get(handles.ms2image,'Value');
if ms2Val == 1
    hAx = handles.axes1;
    %imshow(imadjust(mat2gray(im,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);
    imagesc(im,[0 255] );axis('equal');axis('off'); colormap(gray(256));
end

%meta data:
meta4 = reader.getMetadataStore();
%num of channels (usu 4 for some reason)
nC = meta4.getPixelsSizeC(0).getValue();
%total num of slices:
nSlices = reader.getImageCount();
%number of actual images:
nImages = nSlices/nC;

%set the slider values:
set(handles.slider,'Max',nImages);
set(handles.slider,'Min',1);
%the slider step must be a weird fraction of the max value to ensure integer values:
set(handles.slider,'SliderStep',[1/(nImages-1) 10/(nImages-1)]);
%store fn & fp:
set(handles.ms2dir,'UserData',{fp fn}); 
%report fn to user:
set(handles.ms2dir,'String',fn);
%where to look for the BF reader:
set(handles.load,'UserData',reader);
%set the frame number to 1, for some reason
set(handles.frameNum,'String',1);
set(handles.instrucs,'String','Now use ferature picker to find nuclei or spots. Or load these. Often you have to change the frame to see institute changes');

%gain settings:
mx = 255; %make this easy to change in the future
set(handles.maxSlider,'Min',0);
set(handles.maxSlider,'Max',mx);
set(handles.maxSlider,'Value',mx);
set(handles.maxSlider,'SliderStep',[1/(mx) 10/(mx)]);

set(handles.minSlider,'Min',0);
set(handles.minSlider,'Max',mx);
set(handles.minSlider,'Value',0);
set(handles.minSlider,'SliderStep',[1/(mx) 10/(mx)]);

set(handles.maxScale,'String',mx);
set(handles.minScale,'String',0);

cd(currDir);



function frameNum_Callback(hObject, eventdata, handles)
% hObject    handle to frameNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frameNum as text
%        str2double(get(hObject,'String')) returns contents of frameNum as a double


% --- Executes during object creation, after setting all properties.
function frameNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in select.
function select_Callback(hObject, eventdata, handles)
% hObject    handle to select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%instructions:
% set(handles.instrucs,'String','click and drag to select zoom region');
% 
% aoi = imrect;
% set(handles.instrucs,'String','double click box to confirm');
% %wait for usr to dbl clk:
% wait(aoi);
% set(handles.instrucs,'String','zoom region stored')
% 
% val = round(getPosition(aoi));
% set(handles.aoiValue,'Value',val); %[XMIN YMIN WIDTH HEIGHT]
% 
% hAx = handles.axes1;
% rectangle('Position',val); %,'Parent',hAx);
% 
% str = num2str(val);
% set(handles.aoiValue,'String',str);

%2 pts via mouse clicks
[xpts ypts]=ginput(2);
%get xs & ys, place in order
xpts=round(xpts);ypts=round(ypts);
xpts=sort(xpts);ypts=sort(ypts);
val=[xpts(1) xpts(2) ypts(1) ypts(2)]; 
%input zoom region:
set(handles.aoiValue,'Value',val);
%update display:
magrange_string=['[' num2str(xpts(1)) ' ' num2str(xpts(2)) ' ' num2str(ypts(1)) ' ' num2str(ypts(2)) ']'];
set(handles.aoiValue,'String',magrange_string)



% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
closereq


% --- Executes on button press in ms2image.
function ms2image_Callback(hObject, eventdata, handles)
% hObject    handle to ms2image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ms2image

hisVal = get(handles.hisImage,'Value');
ms2val = get(handles.ms2image,'Value');
if ms2val == 1 
    set(handles.ms2image,'Value',0);
end

% --- Executes on button press in hisImage.
function hisImage_Callback(hObject, eventdata, handles)
% hObject    handle to hisImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hisImage

hisVal = get(handles.hisImage,'Value');
ms2val = get(handles.ms2image,'Value');
if hisVal == 1
    set(handles.ms2image,'Value',0);
end


% --- Executes on button press in load2.
function load2_Callback(hObject, eventdata, handles)
% hObject    handle to load2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%in case the slider has been fucked with:
set(handles.slider,'Value',1);

currDir = cd;
%if we know where the files might be:
try
    dr = get(handles.ms2dir,'UserData');
    cd(dr{1});
end

%get the image file
[fn fp] = uigetfile('*.tiff','pick your histone max proj tif');
reader = bfGetReader([fp fn]);
im = bfGetPlane(reader, 1); %1st image

%plot 1st image to confirm load
hisVal = get(handles.hisImage,'Value');
if hisVal == 1
    hAx = handles.axes1;
    %imshow(imadjust(mat2gray(im,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);
    imagesc(im,[0 255] );axis('equal'); axis('off');colormap(gray(256));
end

%meta data:
meta4 = reader.getMetadataStore();
%num of channels (usu 4 for some reason)
nC = meta4.getPixelsSizeC(0).getValue();
%total num of slices:
nSlices = reader.getImageCount();
%number of actual images:
nImages = nSlices/nC;

%set the slider values:
set(handles.slider,'Max',nImages);
set(handles.slider,'Min',1);
%the slider step must be a weird fraction of the max value to ensure integer values:
set(handles.slider,'SliderStep',[1/(nImages-1) 10/(nImages-1)]);

set(handles.hisDir,'UserData',{fp fn});
%report fn to user:
set(handles.hisDir,'String',fn);
set(handles.load2,'UserData',reader);
set(handles.frameNum,'String',1);
set(handles.instrucs,'String','Now use ferature picker to find nuclei or spots. Or load these. Often you have to change the frame to see institute changes');

%gain settings:
mx = 255; %change this in coordination with mx from load call back
set(handles.maxSlider,'Min',0);
set(handles.maxSlider,'Max',mx);
set(handles.maxSlider,'Value',mx);
set(handles.maxSlider,'SliderStep',[1/(mx) 10/(mx)]);

set(handles.minSlider,'Min',0);
set(handles.minSlider,'Max',mx);
set(handles.minSlider,'Value',0);
set(handles.minSlider,'SliderStep',[1/(mx) 10/(mx)]);

set(handles.maxScale,'String',mx);
set(handles.minScale,'String',0);

cd(currDir);


% --- Executes on button press in crop.
function crop_Callback(hObject, eventdata, handles)
% hObject    handle to crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

aoi = get(handles.aoiValue,'Value'); %[XMIN XMAX YMIN YMAX]


%make sure aoi is specified:
if length(aoi) < 4
    set(handles.instrucs,'Srting','specify an AOI');
else
    %exit MS2 croopped im:
    try
        %get the file reader
        reader = get(handles.load,'UserData');
        omeMeta = reader.getMetadataStore();
        %get the file path:
        dr = get(handles.ms2dir,'UserData');
        fp = dr{1};
        fn = dr{2};
        %write a new file name:
        cropfn = char(extractBefore(fn,'.'));
        fileOutPath = [fp cropfn 'crop.tiff'];
        %number of frames (for some reason they are thought of as z-stacks
        %here:
        nT = omeMeta.getPixelsSizeZ(0).getValue();
        %allocate space
        I = ones(aoi(4),aoi(3),nT); %x & Y flip flopped like in the ome file??
        %pull out the images, crop them:
        for i = 1:nT
            iPlane = reader.getIndex(i - 1, 0, 0) + 1;
            I(:,:,i) = bfGetPlane(reader,iPlane,aoi(1),aoi(2),aoi(3),aoi(4)); %see help bfGetPlane for notation details
        end
        bfsave(I,fileOutPath, 'BigTiff', true); %need big tiff! bfsave is
        ms2flag = 1;
    catch
        ms2flag = 0;
    end
    %exit his cropped im:
    try
        %get the file reader
        Hreader = get(handles.load2,'UserData');
        HomeMeta = Hreader.getMetadataStore();
        %get the file path:
        Hdr = get(handles.hisDir,'UserData');
        Hfp = Hdr{1};
        Hfn = Hdr{2};
        %write a new file name:
        Hcropfn = char(extractBefore(Hfn,'.'));
        HfileOutPath = [Hfp Hcropfn 'crop.tiff'];
        %number of frames (for some reason they are thought of as z-stacks
        %here:
        nT = HomeMeta.getPixelsSizeZ(0).getValue();
        %allocate space
        I = ones(aoi(4),aoi(3),nT); %x & Y flip flopped like in the ome file??
        %pull out the images, crop them:
        for i = 1:nT
            iPlane = Hreader.getIndex(i - 1, 0, 0) + 1;
            I(:,:,i) = bfGetPlane(Hreader,iPlane,aoi(1),aoi(2),aoi(3),aoi(4)); %see help bfGetPlane for notation details
        end
        bfsave(I,HfileOutPath, 'BigTiff', true); %need big tiff! bfsave is
        hisFlag = 1;
    catch
        hisFlag = 0;
    end
end

%update instructions:
if hisFlag == 1 && ms2flag == 1
    set(handles.instrucs,'String','Successfully saved HIS & MS2 cropped images');
elseif hisFlag == 1 && ms2flag == 0
    set(handles.instrucs,'String','saved only HIS cropped image. No MS2 im loaded');
elseif hisFlag == 0 && ms2flag == 1
    set(handles.instrucs,'String','saved only MS2 cropped image. No HIS im loaded');
else
    set(handles.instrucs,'String','no images loaded. nothing saved');
end
        

% --- Executes on button press in showCrop.
function showCrop_Callback(hObject, eventdata, handles)
% hObject    handle to showCrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showCrop


% --- Executes on slider movement.
function minSlider_Callback(hObject, eventdata, handles)
% hObject    handle to minSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
mnVal = get(handles.minSlider,'Value');

set(handles.minScale,'String',mnVal);
% Now update the display
clowval=round(get(handles.minSlider,'Value'));  % set minimum display intensity
chival=round(get(handles.maxSlider,'Value'));   % set maximum display intensity
axes(handles.axes1);
caxis([clowval chival]);    

% --- Executes during object creation, after setting all properties.
function minSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function maxSlider_Callback(hObject, eventdata, handles)
% hObject    handle to maxSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
mxVal = get(handles.maxSlider,'Value');

set(handles.maxScale,'String',mxVal);
% Now update the display
clowval=round(get(handles.minSlider,'Value'));  % set minimum display intensity
chival=round(get(handles.maxSlider,'Value'));   % set maximum display intensity
axes(handles.axes1);
caxis([clowval chival]);    


% --- Executes during object creation, after setting all properties.
function maxSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function maxScale_Callback(hObject, eventdata, handles)
% hObject    handle to maxScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxScale as text
%        str2double(get(hObject,'String')) returns contents of maxScale as a double

 mxx=str2double(get(handles.maxScale,'String'));   
 curr = get(handles.maxSlider,'Value');
%if get(handles.MaxIntensity,'Max') < mxx
                                    % Use the manual input in order to set
                                    % the slider switch maximum for both
                                    % the Min and Max slider swiches
if curr > mxx
    set(handles.maxSlider,'Value',mxx);
end
                                    
set(handles.maxSlider,'Max',mxx); 
set(handles.minSlider,'Max',mxx);

set(handles.maxSlider,'SliderStep',[1/(mxx) 10/(mxx)]);
set(handles.minSlider,'SliderStep',[1/(mxx) 10/(mxx)]);
    
    


% --- Executes during object creation, after setting all properties.
function maxScale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minScale_Callback(hObject, eventdata, handles)
% hObject    handle to minScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minScale as text
%        str2double(get(hObject,'String')) returns contents of minScale as a double


% --- Executes during object creation, after setting all properties.
function minScale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PixelNumber_Callback(hObject, eventdata, handles)
% hObject    handle to PixelNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PixelNumber as text
%        str2double(get(hObject,'String')) returns contents of PixelNumber as a double


% --- Executes during object creation, after setting all properties.
function PixelNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PixelNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FrameRange_Callback(hObject, eventdata, handles)
% hObject    handle to FrameRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrameRange as text
%        str2double(get(hObject,'String')) returns contents of FrameRange as a double


% --- Executes during object creation, after setting all properties.
function FrameRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input_Callback(hObject, eventdata, handles)
% hObject    handle to input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input as text
%        str2double(get(hObject,'String')) returns contents of input as a double



% --- Executes during object creation, after setting all properties.
function input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function output_Callback(hObject, eventdata, handles)
% hObject    handle to output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output as text
%        str2double(get(hObject,'String')) returns contents of output as a double


% --- Executes during object creation, after setting all properties.
function output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in brightnessMinus10.
function brightnessMinus10_Callback(hObject, eventdata, handles)
% hObject    handle to brightnessMinus10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = str2double(get(handles.brightness,'String'));
newVal = num2str(val - 5);
set(handles.brightness,'String',newVal);
%re-pick:
pick_Callback(hObject, eventdata,handles);

% --- Executes on button press in brightnessPlus10.
function brightnessPlus10_Callback(hObject, eventdata, handles)
% hObject    handle to brightnessPlus10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = str2double(get(handles.brightness,'String'));
newVal = num2str(val + 5);
set(handles.brightness,'String',newVal);
%re-pick:
pick_Callback(hObject, eventdata,handles);

% --- Executes on button press in pick.
function pick_Callback(hObject, eventdata, handles)
% hObject    handle to pick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get inputsL
nd = str2double(get(handles.noise,'String'));
b = str2double(get(handles.brightness,'String'));
pixnum = str2double(get(handles.PixelNumber,'String'));
frameNum = str2double(get(handles.frameNum,'String'));
ave = 1; %use this later for frame averageing
shoNum = get(handles.aoiNum,'Value');

%get the image:
im = get(handles.slider,'UserData');

%spot diam must be odd:
sd = round(str2double(get(handles.diameter,'String'))); 
if ~mod(sd,2)
    sd = sd+1;
end

set(handles.pick,'UserData',[]);
%clear the current AOIs:
slider_Callback( handles.frameNum, eventdata, handles) %or something

%im dimensions:
[frmrose frmcol]=size(im);
xlow=1;xhigh=frmcol;ylow=1;yhigh=frmrose; 

%if zoom, restrict limits (come back to this):
crp = get(handles.showCrop,'Value');
if crp == 1
%     limitsxy=eval( get(handles.MagRangeYX,'String') );  % LJF Get the limits of the magnified region
%                                                    % [xlow xhi ylow yhi]
    limitsxy = get(handles.aoiValue,'Value');  %TH
    xlow=limitsxy(1);xhigh=limitsxy(2);            % Define frame limits as those of 
    ylow=limitsxy(3);yhigh=limitsxy(4);            % the magnified region
end

%find spots:
dat=bpass(double(im(ylow:yhigh,xlow:xhigh)),nd,sd);
pk=pkfnd(dat,b,sd);
pk=cntrd(dat,pk,sd+2);

[aoirose aoicol]=size(pk);
                    % Put the aois into our handles structure handles.FitData = [frm#  ave  x   y  pixnum  aoinum]
if aoirose~=0       % If there are spots, put them into handles.FitData and draw them
    pk(:,1)=pk(:,1)+xlow-1;             % Correct coordinates for case where we used a magnified region
    pk(:,2)=pk(:,2)+ylow-1;
    
    FitData=[frameNum*ones(aoirose,1) ave*ones(aoirose,1) pk(:,1) pk(:,2) (pixnum)*ones(aoirose,1) [1:aoirose]'];
                    % Draw the aois
    for indx=1:aoirose
        viscircles(FitData(indx,3:4),pixnum/2,'Color','b','LineWidth',0.5,'EnhanceVisibility',false);
       if shoNum == 1
           text(FitData(indx,3),FitData(indx,4),num2str(FitData(indx,6)),'FontSize',8,'Color','b');
       end
    end
end

if exist('FitData','var')
    set(handles.pick,'UserData',FitData); 
else 
    set(handles.pick,'UserData',[]); %if there are no spots, make this an empty vector.
end

                        %draw_aois(handles.FitData,imagenum,pixnum,handles.DriftList);

function brightness_Callback(hObject, eventdata, handles)
% hObject    handle to brightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of brightness as text
%        str2double(get(hObject,'String')) returns contents of brightness as a double


% --- Executes during object creation, after setting all properties.
function brightness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to brightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in brightnessMinus.
function brightnessMinus_Callback(hObject, eventdata, handles)
% hObject    handle to brightnessMinus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = str2double(get(handles.brightness,'String'));
newVal = num2str(val - 1);
set(handles.brightness,'String',newVal);
%re-pick:
pick_Callback(hObject, eventdata,handles);

% --- Executes on button press in brightnessPlus.
function brightnessPlus_Callback(hObject, eventdata, handles)
% hObject    handle to brightnessPlus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = str2double(get(handles.brightness,'String'));
newVal = num2str(val + 1);
set(handles.brightness,'String',newVal);
%re-pick:
pick_Callback(hObject, eventdata,handles);

function diameter_Callback(hObject, eventdata, handles)
% hObject    handle to diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diameter as text
%        str2double(get(hObject,'String')) returns contents of diameter as a double


% --- Executes during object creation, after setting all properties.
function diameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in diameterMinus.
function diameterMinus_Callback(hObject, eventdata, handles)
% hObject    handle to diameterMinus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = str2double(get(handles.diameter,'String'));
newVal = num2str(val - 2);
set(handles.diameter,'String',newVal);
%re-pick:
pick_Callback(hObject, eventdata,handles);

% --- Executes on button press in diameterPlus.
function diameterPlus_Callback(hObject, eventdata, handles)
% hObject    handle to diameterPlus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = str2double(get(handles.diameter,'String'));
newVal = num2str(val + 2);
set(handles.diameter,'String',newVal);
%re-pick:
pick_Callback(hObject, eventdata,handles);

function noise_Callback(hObject, eventdata, handles)
% hObject    handle to noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noise as text
%        str2double(get(hObject,'String')) returns contents of noise as a double


% --- Executes during object creation, after setting all properties.
function noise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in noiseMinus.
function noiseMinus_Callback(hObject, eventdata, handles)
% hObject    handle to noiseMinus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = str2double(get(handles.noise,'String'));
newVal = num2str(val - 0.1);
set(handles.noise,'String',newVal);
%re-pick:
pick_Callback(hObject, eventdata,handles);

% --- Executes on button press in noisePlus.
function noisePlus_Callback(hObject, eventdata, handles)
% hObject    handle to noisePlus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = str2double(get(handles.noise,'String'));
newVal = num2str(val + 0.1);
set(handles.noise,'String',newVal);
%re-pick:
pick_Callback(hObject, eventdata,handles);

% --- Executes on button press in addAOIs.
function addAOIs_Callback(hObject, eventdata, handles)
% hObject    handle to addAOIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%tell 'em how its done:
set(handles.instrucs,'String','left click spot to add AOI, right click to save and exit')

pixnum = str2double(get(handles.PixelNumber,'String'));
frameNum = str2double(get(handles.frameNum,'String'));
ave = 1; %use this later for frame averageing
%get the image:
% im = get(handles.slider,'UserData');
axes(handles.axes1) 
flag = 0;
aoiinfo = [];
%get the aois, if none picked, set instrucs:
FitData = get(handles.pick,'UserData'); %this may be an empty struc. 
%pick the aois
if ~isempty(FitData)  %if we have spots
    aoinumber=1+max(FitData(:,6));
    while flag == 0
        [a b but]=ginput(1);
        if but==3
            flag=1;
        else
            aoiinfo=[aoiinfo; frameNum ave a b pixnum aoinumber];
            aoinumber=aoinumber+1;                      %Give each aoi a number
            axes(handles.axes1);
            hold on
            [maoi naoi]=size(aoiinfo);
            for indx=maoi:maoi 
                  viscircles(aoiinfo(indx,3:4),pixnum/2,'Color','b','LineWidth',0.5,'EnhanceVisibility',false);
            end
            hold off
        end
    end
    set(handles.pick,'UserData',[FitData; aoiinfo])    
else
    aoinumber = 1;
    while flag == 0
        [a b but]=ginput(1);
        if but==3
            flag=1;
        else
            aoiinfo=[aoiinfo; frameNum ave a b pixnum aoinumber];
            aoinumber=aoinumber+1;                      %Give each aoi a number
            axes(handles.axes1);
            hold on
            [maoi naoi]=size(aoiinfo);
            for indx=maoi:maoi 
                  viscircles(aoiinfo(indx,3:4),pixnum/2,'Color','b','LineWidth',0.5,'EnhanceVisibility',false);
            end
            hold off
        end
    end
    set(handles.pick,'UserData',aoiinfo);
end

set(handles.instrucs,'String','aois added')

% --- Executes on button press in removeAOIs.
function removeAOIs_Callback(hObject, eventdata, handles)
% hObject    handle to removeAOIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%tell 'em how its done:
set(handles.instrucs,'String','left click near spot to remove AOI, right click to save and exit')
%set up the loop:
axes(handles.axes1) 
flag = 0;
%get the aois, if none picked, set instrucs:
FitData = get(handles.pick,'UserData'); %this may be an empty struc. 
while flag==0
    [a b but]=ginput(1);
    if but==3
        flag=1; %right clicking exits the loop
    else
        % Get the aoi number for the aoi closest to where user clicked
        num_closest=aoicompareFli([a b],FitData);
        % logical array, =1 when it matches the aoi number
        logik=(FitData(:,6)==num_closest); 
        % remove information for that aoi
        FitData(logik,:)=[];     
        % Update the existin list of aoi information so that no aoi numbers are skipped:
        FitData=update_FitData_aoinum(FitData);
    end
    % Update the list of AOIs
    set(handles.pick,'UserData',FitData);
    % Refresh display
    slider_Callback(handles.frameNum, eventdata, handles)
end
set(handles.instrucs,'String','AOI list updated')

% --- Executes on button press in findAllSpots.
function findAllSpots_Callback(hObject, eventdata, handles)
% hObject    handle to findAllSpots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Tell user what's up; set up a timer to display at command line
tic;
set(handles.findAllSpots,'String','...')
set(handles.instrucs,'String','finding spots in many frames')
pause(0.1)
%get/set the fp & fn:
%fn & fp from hidden element hisDir or ms2dir:
hisImage = get(handles.hisImage,'Value');
if hisImage == 1
    fileCell = get(handles.hisDir,'UserData');
    fp = fileCell{1}; fn = fileCell{2};
    folder = [fp fn];
else
    fileCell = get(handles.ms2dir,'UserData');
    fp = fileCell{1}; fn = fileCell{2};
    folder = [fp fn];
end
handles.TiffFolder = folder;
guidata(gcbo,handles) ;
%find the spots, now with fitting 
%first check if we do one or two conditions:
cndn = get(handles.useTwo,'Value');
if cndn == 0
    AllSpots=FindAllSpotsFli(handles,3500);  
    set(handles.findAllSpots,'UserData',AllSpots);
elseif cndn == 1
    %start by saving the feature picker inputs to reset later (you'll need
    %to change these as inputs to FindAllSpotsFli):
    noise = get(handles.noise,'String');
    sd = get(handles.diameter,'String');
    sb = get(handles.brightness,'String');
    %find spots under the first condition:
    %to do so we'll set feature picker conditions to what's in the
    %useTwoPanel (bc FindAllSpotsFli uses handles as an input):
    %get 1st cndn inputs:
    noise1 = get(handles.noiseOne,'String');
    sd1 = get(handles.diameterOne,'String');
    sb1 = get(handles.brightnessOne,'String');
    %now reset the feature picker inputs:
    set(handles.noise,'String',noise1)
    set(handles.diameter,'String',sd1)
    set(handles.brightness,'String',sb1)
    %now find spots for cndn 1:
    AllSpots1=FindAllSpotsFli(handles,3500); 
    %repeat for cndn 2:
    noise2 = get(handles.noiseTwo,'String');
    sd2 = get(handles.diameterTwo,'String');
    sb2 = get(handles.brightnessTwo,'String');
    set(handles.noise,'String',noise2)
    set(handles.diameter,'String',sd2)
    set(handles.brightness,'String',sb2)
    AllSpots2=FindAllSpotsFli(handles,3500);
    %compare two allSpots and take the subset of
    %spots that appear in both:
    AllSpots = allSpotsCompare(AllSpots1,AllSpots2);
    %save in gui data:
    set(handles.findAllSpots,'UserData',AllSpots);
    %return feature picker settings to previous values:
    set(handles.noise,'String',noise);
    set(handles.diameter,'String',sd);
    set(handles.brightness,'String',sb);
end
toc
%update:
set(handles.instrucs,'String','spots found. now you can save them.')
set(handles.findAllSpots,'String','Find all spots')

% --- Executes on button press in saveAllSpots.
function saveAllSpots_Callback(hObject, eventdata, handles)
% hObject    handle to saveAllSpots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fl = get(handles.outputText,'UserData'); % file locations are loaded at GUI launch
filestring = get(handles.output,'String');
AllSpots = get(handles.findAllSpots,'UserData');

eval(['save ' fl.data filestring ' AllSpots']);

set(handles.output,'String','default.dat')
set(handles.instrucs,'String','spots saved');

% --- Executes on button press in loadAllSpots.
function loadAllSpots_Callback(hObject, eventdata, handles)
% hObject    handle to loadAllSpots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fl = get(handles.outputText,'UserData'); % file locations are loaded at GUI launch
filestring = get(handles.input,'String');

try
    eval(['load ' fl.data filestring ' -mat']) %this should load the struc 'AllSpots'
catch
    eval(['load ' fl.outPath filestring ' -mat'])
end

set(handles.findAllSpots,'UserData',AllSpots);

set(handles.instrucs,'String','spots loaded');
set(handles.input,'String','default.dat')

% --- Executes on button press in loadNuclei.
function loadNuclei_Callback(hObject, eventdata, handles)
% hObject    handle to loadNuclei (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fl = get(handles.outputText,'UserData'); % file locations are loaded at GUI launch
filestring = get(handles.input,'String');

eval(['load ' fl.data filestring ' -mat'])

%get what type:
if exist('aoiinfo2','var')
    set(handles.pick,'UserData',aoiinfo2);
    slider_Callback( handles.frameNum, eventdata, handles)
    set(handles.instrucs,'String','loaded nuclei locations. not tracking data.') 
elseif exist('aoifits','var')
    set(handles.pick,'UserData',aoifits.aoiinfo2);
    slider_Callback( handles.frameNum, eventdata, handles)
    set(handles.instrucs,'String','loaded nuclei location & tracking data.')
    set(handles.trackNuclei,'UserData',aoifits);
elseif exist('AllSpots','var')
    set(handles.pick,'UserData',AllSpots.aoiinfo2);
    slider_Callback( handles.frameNum, eventdata, handles)
    %need to do someting with the rest of AllSpots
else
    set(handles.instrucs,'String','unrecognized file type')    
end

set(handles.input,'String','default.dat')
   
% --- Executes on button press in saveNuclei.
function saveNuclei_Callback(hObject, eventdata, handles)
% hObject    handle to saveNuclei (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fl = get(handles.outputText,'UserData'); % file locations are loaded at GUI launch
filestring = get(handles.output,'String');
aoiinfo2 = get(handles.pick,'UserData');
nukeNum = aoiinfo2(end,end);
fprintf('%d nuclei saved \n',nukeNum);

eval(['save ' fl.data filestring ' aoiinfo2']);

set(handles.output,'String','default.dat')
update = sprintf('%d nuclei locations saved',nukeNum);
set(handles.instrucs,'String',update);
%plot in yellow:
aoirose = size(aoiinfo2,1);
pixnum = str2double(get(handles.PixelNumber,'String'));
for indx=1:aoirose
    viscircles(aoiinfo2(indx,3:4),pixnum/2,'Color','y','LineWidth',0.5,'EnhanceVisibility',false);
end  

% --- Executes on button press in trackNuclei.
function trackNuclei_Callback(hObject, eventdata, handles)
% hObject    handle to trackNuclei (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tic;
%set(handles.trackNuclei,'String','...')
set(handles.instrucs,'String','tracking. this may take a moment. 0.0017s /frame /pix /nuclei on Hardens machine')
pause(0.1)

%some inputs:
pixnum = str2double(get(handles.PixelNumber,'String'));
frms = eval([get(handles.FrameRange,'String') ';']);
ave = 1; %use this later for frame averageing
aoiinf = get(handles.pick,'UserData');
FitData = get(handles.pick,'UserData');
%for tracking forwards and backwards
pickFrame = round(mean(aoiinf(:,1))); %fine the frame when (most) of the aois were selected in

%fn & fp from hidden element hisDir or ms2dir:
hisImage = get(handles.hisImage,'Value');
if hisImage == 1
    fileCell = get(handles.hisDir,'UserData');
    fp = fileCell{1}; fn = fileCell{2};
    folder = [fp fn];
else
    fileCell = get(handles.ms2dir,'UserData');
    fp = fileCell{1}; fn = fileCell{2};
    folder = [fp fn];
end
handles.TiffFolder = folder;
guidata(gcbo,handles) ; %need to update handles for DS fns gauss2d+mapstruc2d_v2Fli and friends

%from LJF
[mfrms nfrms]=size(frms);
if mfrms ==1
    frms=frms';                                 % frms now a column vector
end
% %to reverse track: %191219 we unimplemented reverse tracking
% rev = get(handles.reverse,'Value');
% if rev == 1
%     frms = flip(frms);
% end

%aoi stuff:
[maoi naoi]=size(aoiinf); % maoi is the number of aois
%start fitting stuff (LJF):
argoutsImageData=[];
argoutsBackgroundData=[];

aoifits.dataDescription='[aoinumber framenumber amplitude xcenter ycenter sigma offset integrated_aoi (integrated pixnum) (original aoi#)]';
aoifits.parameter=[ave pixnum];
aoifits.parameterDescription='[ave aoiDiameter]';
aoifits.tifFile=folder;
aoifits.centers=aoiinf(:,3:4);
aoifits.centersDescription='[aoi_xcenters aoi_ycenters]';
aoifits.aoiinfo2Description='[(framenumber when marked) ave x y aoiDiameter aoinumber]';
aoifits.aoiinfo2=aoiinf;  %may need to inspect/change
aoifits.AllSpotsDescription='aoifits.AllSpots{m,1}=[x y] spots in frm m;  {m,2}=# of spots,frame m; {m,3}=frame #; {1,4}=[firstframe:lastframe]; {2,4}=NoiseDiameter  SpotDiameter  SpotBrightness]'  ;
                                        '{2,4}=NoiseDiameter  SpotDiameter  SpotBrightness]'  ;
% aoifits.AllSpots=FreeAllSpotsMemory(handles.AllSpots);      %We need to make an AllSpost struc, then can plug it into here. we've yet to do this
%outputName=get(handles.output,'String'); %presently we will do this in saveTracks button

%determine the nature of tracking
frmIndx = find(frms == pickFrame); %gives the index of the frame we picked the nukes in rel. to the frame range we track over

% Build a 2D mapstruc to direct data processing, also check for existance
% of DL:
useDL = get(handles.useDriftlist,'Value');
if useDL == 0  %no DL
    if frmIndx == 1 %forward tracking
        mapstruc2d=build_2d_mapstruc_aois_frmsFli(folder,pixnum,frms,ave,FitData,1,handles); 
        DataOutput2d=gauss2d_mapstruc2d_v2Fli(mapstruc2d,handles); % Process the data (integrate, fit etc) V.2 is parallel processing  
    elseif frmIndx == nfrms %reverse tracking
        frms = flip(frms);
        mapstruc2d=build_2d_mapstruc_aois_frmsFli(folder,pixnum,frms,ave,FitData,1,handles);
        DataOutput2d=gauss2d_mapstruc2d_v2Fli(mapstruc2d,handles); % Process the data (integrate, fit etc) V.2 is parallel processing  
    else  %track forward & back from the pickFrame
        %first do the reverse tracking
        firstFrms = flip(frms(1:frmIndx));
        mapstruc2dFirst=build_2d_mapstruc_aois_frmsFli(folder,pixnum,firstFrms,ave,FitData,1,handles);
        DataOutput2dFirst=gauss2d_mapstruc2d_v2Fli(mapstruc2dFirst,handles); % Process the data (integrate, fit etc) V.2 is parallel processing 
        %%% we need to figure out how to re-order the above two outputs to
        %%% prevent a downstrewam fuckup
        %%%
        %now the forward tracking
        lastFrms = frms(frmIndx + 1:nfrms);
        mapstruc2dLast=build_2d_mapstruc_aois_frmsFli(folder,pixnum,lastFrms,ave,FitData,1,handles);
        DataOutput2dLast=gauss2d_mapstruc2d_v2Fli(mapstruc2dLast,handles); % Process the data (integrate, fit etc) V.2 is parallel processing  
        %combine the two mats in the DataOutput strucs:
        DataOutput2d.ImageData = [DataOutput2dFirst.ImageData; DataOutput2dLast.ImageData];
        DataOutput2d.BackgroundData = [DataOutput2dFirst.BackgroundData; DataOutput2dLast.BackgroundData];
    end        
else
    %for DL, we only do forward tracking for now
    %make sure we have a DL
    DL = get(handles.createDriftlist,'UserData');
    if isempty(DL)
        set(handles.instrucs,'String','No drift list loaded. Make or load that then try again')
        return;
    else
        mapstruc2d=build_2d_mapstruc_aois_frmsFli(folder,pixnum,frms,ave,FitData,2,handles);        % use a measure of global nuclei drift (a DL) to track aois. its mo betta
        %all the DL magic in the tracking pipeline takes place in the
        %sub-script build_mapstruc_cell_columnFli
        DataOutput2d=gauss2d_mapstruc2d_v2Fli(mapstruc2d,handles); % Process the data (integrate, fit etc) V.2 is parallel processing  
    end 
end

argoutsImageData=DataOutput2d.ImageData;
argoutsBackgroundData=DataOutput2d.BackgroundData;
     % Start a gui for display of the fit results
      % argouts=[ aoinumber framenumber amplitude xcenter ycenter sigma offset integrated_aoi]
       % Save the data after each aoi is processed  
                % First assign the ImageData

%eval(['save p:\matlab12\larry\data\' outputName ' aoifits' ])
aoifits.data=argoutsImageData;
aoifits.BackgroundData=argoutsBackgroundData;

%%% HERES WHERE WE PUT AOIFITS %%%
set(handles.trackNuclei,'UserData',aoifits);
%reset strings:
set(handles.instrucs,'String','nuclei tracked')
%set(handles.trackNuclei,'String','Track !')
toc

% --- Executes on button press in loadNucleiTracks.
function loadNucleiTracks_Callback(hObject, eventdata, handles)
% hObject    handle to loadNucleiTracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fl = get(handles.outputText,'UserData'); % file locations are loaded at GUI launch
filestring = get(handles.input,'String');

try
    eval(['load ' fl.data filestring ' -mat'])
catch
    eval(['load ' fl.outPath filestring ' -mat'])
end

%get what type:
if exist('aoifits','var')
    set(handles.pick,'UserData',aoifits.aoiinfo2);
    slider_Callback( handles.frameNum, eventdata, handles)
    set(handles.instrucs,'String','loaded nuclei tracking data.')
    set(handles.trackNuclei,'UserData',aoifits);
elseif exist('aoiinfo2','var')
    set(handles.instrucs,'String','no tracking data in this .dat file. only nuclei locations loaded.')
    set(handles.pick,'UserData',aoiinfo2);
    slider_Callback( handles.frameNum, eventdata, handles)
% elseif AllSpots
%     set(handles.pick,'UserData',AllSpots.aoiinfo2);
%     slider_Callback( handles.frameNum, eventdata, handles)
%     %need to do someting with the rest of AllSpots
else
    set(handles.instrucs,'String','unrecognized file type')    
end

set(handles.input,'String','default.dat')

% --- Executes on button press in saveNucleiTracks.
function saveNucleiTracks_Callback(hObject, eventdata, handles)
% hObject    handle to saveNucleiTracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

aoifits = get(handles.trackNuclei,'UserData');
fl = get(handles.outputText,'UserData'); % file locations are loaded at GUI launch
filestring = get(handles.output,'String');

eval(['save ' fl.data filestring ' aoifits']);

set(handles.output,'String','default.dat')
set(handles.instrucs,'String','nuclei tracks saved');

function EditRollingBallRadius_Callback(hObject, eventdata, handles)
% hObject    handle to EditRollingBallRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditRollingBallRadius as text
%        str2double(get(hObject,'String')) returns contents of EditRollingBallRadius as a double


% --- Executes during object creation, after setting all properties.
function EditRollingBallRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditRollingBallRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditRollingBallHeight_Callback(hObject, eventdata, handles)
% hObject    handle to EditRollingBallHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditRollingBallHeight as text
%        str2double(get(hObject,'String')) returns contents of EditRollingBallHeight as a double


% --- Executes during object creation, after setting all properties.
function EditRollingBallHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditRollingBallHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in IncrementRollingBallRadius.
function IncrementRollingBallRadius_Callback(hObject, eventdata, handles)
% hObject    handle to IncrementRollingBallRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%grab curr val:
val = str2double(get(handles.EditRollingBallRadius,'String'));
newVal = num2str(val + 2);
%set new val:
set(handles.EditRollingBallRadius,'String',newVal);
%re-set display:
slider_Callback( handles.frameNum, eventdata, handles)

% --- Executes on button press in DecrementRollingBallRadius.
function DecrementRollingBallRadius_Callback(hObject, eventdata, handles)
% hObject    handle to DecrementRollingBallRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%grab curr val:
val = str2double(get(handles.EditRollingBallRadius,'String'));
newVal = num2str(val - 2);
%set new val:
set(handles.EditRollingBallRadius,'String',newVal);
%re-set display:
slider_Callback( handles.frameNum, eventdata, handles)

% --- Executes on button press in IncrementRollingBallHeight.
function IncrementRollingBallHeight_Callback(hObject, eventdata, handles)
% hObject    handle to IncrementRollingBallHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%grab curr val:
val = str2double(get(handles.EditRollingBallHeight,'String'));
newVal = num2str(val + 1);
%set new val:
set(handles.EditRollingBallHeight,'String',newVal);
%re-set display:
slider_Callback( handles.frameNum, eventdata, handles)

% --- Executes on button press in DecrementRollingBallHeight.
function DecrementRollingBallHeight_Callback(hObject, eventdata, handles)
% hObject    handle to DecrementRollingBallHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%grab curr val:
val = str2double(get(handles.EditRollingBallHeight,'String'));
newVal = num2str(val - 1);
%set new val:
set(handles.EditRollingBallHeight,'String',newVal);
%re-set display:
slider_Callback( handles.frameNum, eventdata, handles)

% --- Executes on selection change in BackgroundChoice.
function BackgroundChoice_Callback(hObject, eventdata, handles)
% hObject    handle to BackgroundChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BackgroundChoice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BackgroundChoice


% --- Executes during object creation, after setting all properties.
function BackgroundChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BackgroundChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in FitChoice.
function FitChoice_Callback(hObject, eventdata, handles)
% hObject    handle to FitChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FitChoice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FitChoice


% --- Executes during object creation, after setting all properties.
function FitChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FitChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadMatch.
function loadMatch_Callback(hObject, eventdata, handles)
% hObject    handle to loadMatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fl = get(handles.outputText,'UserData'); % file locations are loaded at GUI launch
filestring = get(handles.input,'String');

try
    eval(['load ' fl.data filestring ' -mat'])
catch
    eval(['load ' fl.outPath filestring ' -mat'])
end

%what are we calling the match output? ...store it in makeMatches user data
set(handles.makeMatches,'UserData',analData);
% get the tracks info: ...incomplete. make changes to lines 255 in slider callback
% aoifits = analData.nukesMat;
% set(handles.trackNuclei,'UserData',aoifits);

%old:
% %get what type:
% if exist('aoifits','var')
%     set(handles.pick,'UserData',aoifits.aoiinfo2);
%     slider_Callback( handles.frameNum, eventdata, handles)
%     set(handles.instrucs,'String','loaded nuclei locations. not tracking data.')
%     set(handles.trackNuclei,'UserData',aoifits);
% elseif exist('aoiinfo2','var')
%     set(handles.instrucs,'String','not tracking data here. only locations loaded.')
%     set(handles.pick,'UserData',aoiinfo2);
%     slider_Callback( handles.frameNum, eventdata, handles)
% % elseif AllSpots
% %     set(handles.pick,'UserData',AllSpots.aoiinfo2);
% %     slider_Callback( handles.frameNum, eventdata, handles)
% %     %need to do someting with the rest of AllSpots
% else
%     set(handles.instrucs,'String','unrecognized file type')    
% end

set(handles.input,'String','default.dat')
set(handles.instrucs,'String','paired spots loaded')


% --- Executes on button press in makeMatches.
function makeMatches_Callback(hObject, eventdata, handles)
% hObject    handle to makeMatches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tic;
%get tracked nuclei aoifits struc
nukes = get(handles.trackNuclei,'UserData');
%get spots AllSpots struc
spots = get(handles.findAllSpots,'UserData');
%time base vector:
timebaseStruc = get(handles.timeBase,'UserData');
timeBase = timebaseStruc.ttb;
if isempty(timeBase)
    set(handles.instrucs,'String','No time base. Make one and try again.')
    set(handles.makeMatches,'String','Pair spots')
    return
end
%get nuclear cycle info:
ncs = get(handles.nuclearCycles,'UserData');
if isempty(ncs)
    set(handles.instrucs,'String','No NC input. Add that and try again.')
    set(handles.makeMatches,'String','Pair spots')
    return
end
%call match making script:
if isempty(timeBase)
    set(handles.instrucs,'String','No time base, continuing without it')
    set(handles.makeMatches,'String','...')
    pause(0.1)
    analData = matchMaker(nukes,spots,timeBase,ncs);
else 
    set(handles.instrucs,'String','Pairing spots with nuclei')
    set(handles.makeMatches,'String','...')
    pause(0.1)  
    analData = matchMaker(nukes,spots,timeBase,ncs);
end

set(handles.makeMatches,'UserData',analData)
set(handles.makeMatches,'String','Pair spots')
set(handles.instrucs,'String','Spots paired to nuclei. analData is now ready to be saved and analyzed. AP pos. addition is optional.')
toc

%       nuclei locations: handles.pick,'UserData' (aoiinfo2 mat)  
%       nuclei tracking: handles.trackNuclei,'UserData' (aoifits struc)
%       spots: handles.findAllSpots,'UserData' (AllSpots struc)
%       spot-nuclei pairs: handles.makeMatches,'UserData' (this is presently a pairedSpots mat, but this will change)


% --- Executes on button press in saveMatch.
function saveMatch_Callback(hObject, eventdata, handles)
% hObject    handle to saveMatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

analData = get(handles.makeMatches,'UserData');
fl = get(handles.outputText,'UserData'); % file locations are loaded at GUI launch
filestring = get(handles.output,'String');

eval(['save ' fl.data filestring ' analData']);

set(handles.output,'String','default.dat')
set(handles.instrucs,'String','analData structure saved');

% --- Executes on button press in plotSpots.
function plotSpots_Callback(hObject, eventdata, handles)
% hObject    handle to plotSpots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotSpots


% --- Executes on button press in aoiNum.
function aoiNum_Callback(hObject, eventdata, handles)
% hObject    handle to aoiNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of aoiNum


% --- Executes on button press in removeClose.
function removeClose_Callback(hObject, eventdata, handles)
% hObject    handle to removeClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get the radius and current set of nuclei:
rad = str2double(get(handles.editRadius,'String'));
FitData = get(handles.pick,'UserData');
%get rid of close ones, reset the aoi list:
newFitData = Remove_Close_AOIs_v1Fli(FitData,rad);
set(handles.pick,'UserData',newFitData)
%plot new nuclei (first turn on display features)
set(handles.displayNuclei,'Value',1)
slider_Callback( handles.frameNum, eventdata, handles) 
%slider_Callback( handles, eventdata, handles) 

function editRadius_Callback(hObject, eventdata, handles)
% hObject    handle to editRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editRadius as text
%        str2double(get(hObject,'String')) returns contents of editRadius as a double


% --- Executes during object creation, after setting all properties.
function editRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in editRadiusMinus.
function editRadiusMinus_Callback(hObject, eventdata, handles)
% hObject    handle to editRadiusMinus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = str2double(get(handles.editRadius,'String'));
newVal = num2str(val - 1);
set(handles.editRadius,'String',newVal);

% --- Executes on button press in editRadiusPlus.
function editRadiusPlus_Callback(hObject, eventdata, handles)
% hObject    handle to editRadiusPlus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = str2double(get(handles.editRadius,'String'));
newVal = num2str(val + 1);
set(handles.editRadius,'String',newVal);



function frameFreq_Callback(hObject, eventdata, handles)
% hObject    handle to frameFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frameFreq as text
%        str2double(get(hObject,'String')) returns contents of frameFreq as a double


% --- Executes during object creation, after setting all properties.
function frameFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in displayNuclei.
function displayNuclei_Callback(hObject, eventdata, handles)
% hObject    handle to displayNuclei (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of displayNuclei


% --- Executes on button press in timeBase.
function timeBase_Callback(hObject, eventdata, handles)
% hObject    handle to timeBase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%old way, from user input using startTime & frameFreq callbacks (crrently hidden):
% %get starting time:
% startT = str2double(get(handles.startTime,'String'));
% %get frame frequencey:
% freak = str2double(get(handles.frameFreq,'String'));
% %Get number of frames;
% mxFrame = get(handles.slider,'Max');
% %calc the max time;
% mxTime = mxFrame*freak + (startT - freak); %this avoids the fencepost error and I am pretty sure its correct
% %make time base vector:
% timeBase = startT:freak:mxTime; %step by freak from startT to the max time
% set(handles.timeBase,'UserData',timeBase);
% set(handles.timeBase,'String','Remake time base');
% set(handles.instrucs,'String','time base created. will be added to output once spots are paired.');

%new way, mining the metadata;
set(handles.instrucs,'String','making time base. This may take a minute.')
pause(0.1);

timebaseStruc = makeTimeBaseLS;
set(handles.timeBase,'UserData',timebaseStruc);
set(handles.instrucs,'String','time base created. will be added to output once spots are paired.');

function startTime_Callback(hObject, eventdata, handles)
% hObject    handle to startTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startTime as text
%        str2double(get(hObject,'String')) returns contents of startTime as a double


% --- Executes during object creation, after setting all properties.
function startTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in restrictStepSize.
function restrictStepSize_Callback(hObject, eventdata, handles)
% hObject    handle to restrictStepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of restrictStepSize

%this does not noticably slow down tracking step
val = get(handles.restrictStepSize,'Value');
if val == 1
    set(handles.maxStep,'Visible','On')
    set(handles.maxStepText,'Visible','On')
else
    set(handles.maxStep,'Visible','Off')
    set(handles.maxStepText,'Visible','Off')
end

function maxStep_Callback(hObject, eventdata, handles)
% hObject    handle to maxStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxStep as text
%        str2double(get(hObject,'String')) returns contents of maxStep as a double


% --- Executes during object creation, after setting all properties.
function maxStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in displayTracks.
function displayTracks_Callback(hObject, eventdata, handles)
% hObject    handle to displayTracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of displayTracks


% --- Executes on button press in manualPick.
function manualPick_Callback(hObject, eventdata, handles)
% hObject    handle to manualPick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%initialize:
aoiinfo=[];
frameNum = str2double(get(handles.frameNum,'String'));
ave=1;
pixnum=str2double(get(handles.PixelNumber,'String'));
shoNum = get(handles.aoiNum,'Value');

flag=0;
axes(handles.axes1)
aoinumber=1;
while flag==0
    [a b but]=ginput(1);
    if but==3
        flag=1;
    else
        aoiinfo=[aoiinfo; frameNum ave a b pixnum aoinumber];
        aoinumber=aoinumber+1;                      %Give each aoi a number
        %axes(handles.axes1);
        hold on
        [maoi ~]=size(aoiinfo);
        for indx=maoi:maoi                        
                    % Here to draw aoi boxes only at pixel boundaries
           viscircles(aoiinfo(indx,3:4),pixnum/2,'Color','b','LineWidth',0.5,'EnhanceVisibility',false);
           if shoNum == 1
               text(aoiinfo(indx,3),aoiinfo(indx,4),num2str(aoiinfo(indx,6)),'FontSize',8,'Color','b');
           end
        end
        hold off
    end
end
FitData = aoiinfo; %only do this to let user know that this is the same thing you get from auto picker
set(handles.pick,'UserData',FitData);
% handles.FitData=aoiinfo;                            % Store the list in the handles structure
% guidata(gcbo,handles) ;


% --- Executes on button press in createDriftlist.
function createDriftlist_Callback(hObject, eventdata, handles)
% hObject    handle to createDriftlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%time base (make sure we have one):
timebaseStruc = get(handles.timeBase,'UserData');
vid.ttb = timebaseStruc.ttb; %originally ttb was a row vector in some absolute time (ms) so you have to subtract the value of the first mamber. in flimscroll, its also a col vector, in units of (s), with whatever the user puts as time zero. 
if isempty(vid.ttb) 
    set(handles.instrucs,'String','drift list aborted. no time base. make that then try again.')
    pause(0.1);
    return;
end
%use tracking button to track AOIs, but first make sure the 'Use drift list
%button value is zero
set(handles.useDriftlist,'Value',0)
trackNuclei_Callback(handles, eventdata, handles)
% trackNuclei_Callback(hObject, eventdata, handles)
% set(handles.instrucs,'String','tracking nuclei for DL creation.')
% pause(0.1);
%get the results of the tracking
set(handles.instrucs,'String','nuclei tracked. Now making drift list.')
pause(0.1);
FitData = get(handles.trackNuclei,'UserData');

%now figure out how to make a DL...
%1. create dat:
dat = draw_aoifits_aois_v1Fli(FitData);
%2. inputs for DL construction:
SequenceLength = get(handles.slider,'Max');
frms = eval([get(handles.FrameRange,'String') ';']);
CorrectionRange = [frms(1) frms(end)];
%make xy_cell cell array (this is full of vestiges from LJF stuff, but I didn;t get around to editing construct_driftlist_time_v1 tto thoroughly and dont think it slows things much):
trackedAOInum = size(dat,3);
for i = 1:trackedAOInum
    xy_cell{i}.dat = dat(:,:,i);
    xy_cell{i}.range = CorrectionRange;
    xy_cell{i}.userange = CorrectionRange;
end
%3. make DL:
drifts = construct_driftlist_time_v1Fli(xy_cell,vid,CorrectionRange,SequenceLength,[22 22],[3 5 3 5]); %(xy_cell,vid,CorrectionRange,SequenceLength,Polyorderxy,varargin)
%the [22 22] poly order bullshit is not used, as we opt to use the optional
%SG smoothing. 
% DL = drifts.diffdriftlist;
DL = drifts; % update for LS by TH 220427

set(handles.createDriftlist,'UserData',DL)
set(handles.instrucs,'String','Driftlist created.')

% --- Executes on button press in saveDriftlist.
function saveDriftlist_Callback(hObject, eventdata, handles)
% hObject    handle to saveDriftlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fl = get(handles.outputText,'UserData'); % file locations are loaded at GUI launch
filestring = get(handles.output,'String');
DL = get(handles.createDriftlist,'UserData');

eval(['save ' fl.data filestring ' DL']);

set(handles.output,'String','default.dat')
set(handles.instrucs,'String','Drift list saved');


% --- Executes on button press in loadDriftlist.
function loadDriftlist_Callback(hObject, eventdata, handles)
% hObject    handle to loadDriftlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fl = get(handles.outputText,'UserData'); % file locations are loaded at GUI launch
filestring = get(handles.input,'String');

eval(['load ' fl.data filestring ' -mat'])

set(handles.createDriftlist,'UserData',DL);
set(handles.instrucs,'String','Drift list loaded.')
set(handles.input,'String','default.dat')

% --- Executes on button press in useDriftlist.
function useDriftlist_Callback(hObject, eventdata, handles)
% hObject    handle to useDriftlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useDriftlist


% --- Executes on button press in findAP.
function findAP_Callback(hObject, eventdata, handles)
% hObject    handle to findAP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%see this video to figure out how to get an output from a gui: https://blogs.mathworks.com/videos/2010/02/12/advanced-getting-an-output-from-a-guide-gui/
%also this site: https://www.mathworks.com/matlabcentral/answers/326124-how-to-manage-more-than-one-gui-ie-main-gui-and-sub-gui
%from the latter, we get:

%call painInTheAxis:
coordStruc = painInTheAxis;
%maybe make pain in the axix modal?

if ~isempty(coordStruc)
    %we get the output from pain in the axis and save it
    %here
    set(handles.findAP,'UserData',coordStruc)
    %the APpos should also be saved to a .dat file by painintheaxis for easy
    %recall.
    set(handles.instrucs,'String','A-P position found, loaded and saved.')
else
    set(handles.instrucs,'String','A-P position interrupted. AP not saved.')
end

% --- Executes on button press in loadAP.
function loadAP_Callback(hObject, eventdata, handles)
% hObject    handle to loadAP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currDir = cd;
%if a file has already been retrieved, go to that dir
try
    dr = get(handles.hisDir,'UserData');
    cd(dr{1});
    cd ..
end

%get the image file
[fn fp] = uigetfile('*.dat','pick your AP position dat file');
if ischar(fn)  %load, if uninterrupted  
    eval(['load ' fp fn ' -mat'])
    set(handles.findAP,'UserData',coordStruc)

    set(handles.instrucs,'String','AP position file loaded. Now add AP positions to analData struc.')
    cd(currDir);
else
    set(handles.instrucs,'String','interrupted. no AP location loaded')
    cd(currDir);
    return
end

% --- Executes on button press in addAPpos.
function addAPpos_Callback(hObject, eventdata, handles)
% hObject    handle to addAPpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get AP coords:
coordStruc = get(handles.findAP,'UserData');
try 
    coordA = coordStruc.coordA; %[xpos ypos]
    coordP = coordStruc.coordP;
catch
    set(handles.instrucs,'String','no AP coords saved. get that and come back')
    return
end
%get the spots thingy:
analData = get(handles.makeMatches,'UserData');
try 
    spotsMat = analData.spotsMat; %[1.time 2.frameNum 3.nucleus 4.distanceFromNucleus 5.xPos 6.yPos 7.amplitude 8.sigma 9.offset 10.intensity 11.relativeFrameNumber]';
    nukesMat = analData.nukesMat;
catch
    set(handles.instrucs,'String','no AP coords saved. get that and come back')
    return
end
%geometry:
%angle:
APangle=atan2((-coordP(2)+coordA(2)),(coordP(1)-coordA(1)));    %in radians
%TH use this to avoid a sign error in the angle calc. I cannot imagine how this
%matters, and may make things worse? The problem stems from the fact that
%the y-coords are negative in the ref frame of the embryos. so when you
%select the A or P position, you get a negative y-value.
% length:
APlength=sqrt((-coordP(1)+coordA(1))^2+(-coordP(2)+coordA(2))^2); %in pixels
%convert spot locations:
angles=atan2((-spotsMat(:,6)+coordA(2)),(spotsMat(:,5)-coordA(1)));
%Distance between the points and the ANT point:
distances=sqrt((coordA(2)-spotsMat(:,6)).^2+(coordA(1)-spotsMat(:,5)).^2);
%get dist along AP axis:
APpositions=distances.*cos(angles-APangle);
%add to mat:
spotsMat(:,13)=APpositions./APlength; %check for sign errors in the future. they will likely arise from 'angles'. soln uis likely just to make this 'abs(... '
%repeat for nuclei:
nukeAngles=atan2((-nukesMat(:,6)+coordA(2)),(nukesMat(:,5)-coordA(1)));
%Distance between the points and the ANT point:
nukeDistances=sqrt((coordA(2)-nukesMat(:,6)).^2+(coordA(1)-nukesMat(:,5)).^2);
%get dist along AP axis:
nukeAPpositions=nukeDistances.*cos(nukeAngles-APangle);
%add to mat:
nukesMat(:,13)=nukeAPpositions./APlength; %check for sign errors in the future. they will likely arise from 'angles'. soln uis likely just to make this 'abs(... '

%update descriptions:
analData.spotsMatDescription = '[1.time(s) 2.frameNum 3.nucleus 4.distanceFromNucleusCntr 5.xPos 6.yPos 7.gaussAmp 8.sigma 9.offset 10.integratedIntensity 11.relativeFrameNumber 12.NC 13.APpos]';
analData.nukesMatDescription = '[1.time(s) 2.frameNum 3.nucleus 4.binarySpot 5.nukeXpos 6.nukeYpos 7.spotGaussAmp 8.spotSigma 9.spotOffset 10.spotIntegratedIntensity 11.relativeFrameNumber 12.NC 13.nuclearAPpos]';
%save to GUI data
analData.spotsMat = spotsMat;
analData.nukesMat = nukesMat;
set(handles.makeMatches,'UserData',analData)

set(handles.instrucs,'String','AP positions added. analData is now good to go. can save it and exit.')


% --- Executes on button press in saveAPanal.
function saveAPanal_Callback(hObject, eventdata, handles)
% hObject    handle to saveAPanal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

saveMatch_Callback(hObject, eventdata, handles)

% --- Executes on button press in displayDriftlist.
function displayDriftlist_Callback(hObject, eventdata, handles)
% hObject    handle to displayDriftlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of displayDriftlist


% --- Executes on selection change in nc.
function nc_Callback(hObject, eventdata, handles)
% hObject    handle to nc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns nc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nc


% --- Executes during object creation, after setting all properties.
function nc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ncFrame_Callback(hObject, eventdata, handles)
% hObject    handle to ncFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ncFrame as text
%        str2double(get(hObject,'String')) returns contents of ncFrame as a double


% --- Executes during object creation, after setting all properties.
function ncFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ncFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nuclearCycles.
function nuclearCycles_Callback(hObject, eventdata, handles)
% hObject    handle to nuclearCycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get the frame number:
cycleFrame = str2double(get(handles.ncFrame,'String'));
%see if we've added any previous ncs:
ncs = get(handles.nuclearCycles,'UserData');
if isempty(ncs)
    %Initialize nc mat;
    ncs = [10:14]';
    ncs(:,2) = NaN;
end

switch get(handles.nc,'Value')
    case 1
        ncs(5,2) = cycleFrame;
        set(handles.instrucs,'String','Added start of NC 14')
        currNc = 14;
    case 2
        ncs(4,2) = cycleFrame;
        set(handles.instrucs,'String','Added start of NC 13')
        currNc = 13;
    case 3
        ncs(3,2) = cycleFrame;
        set(handles.instrucs,'String','Added start of NC 12')
        currNc = 12;
    case 4
        ncs(2,2) = cycleFrame;
        set(handles.instrucs,'String','Added start of NC 11')
        currNc = 11;
    case 5
        ncs(1,2) = cycleFrame;
        set(handles.instrucs,'String','Added start of NC 10')
        currNc = 10;
end

set(handles.nuclearCycles,'UserData',ncs)
set(handles.ncFrame,'String','--')
%update at the command line too:
fprintf('Added frame %d as start of NC %d \n',cycleFrame,currNc);


% --- Executes on button press in saveTime.
function saveTime_Callback(hObject, eventdata, handles)
% hObject    handle to saveTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

timebaseStruc = get(handles.timeBase,'UserData');
%jump around directories then save:
currDir = cd;
dr = get(handles.ms2dir,'UserData');
cd(dr{1})
cd ..
fp = cd;
fn = extractBefore(dr{2},'_'); %try this, should work well enough to get some info on the image to the filename    
fileOutPath = [fp filesep fn '-ttb.dat'];
eval(['save ' fileOutPath ' timebaseStruc']);
cd(currDir)

set(handles.instrucs,'String','Saved time base.')
fprintf('saved time base to %s',fileOutPath)

% --- Executes on button press in loadTime.
function loadTime_Callback(hObject, eventdata, handles)
% hObject    handle to loadTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%navigate to where the file likely will be:
currDir = cd;
dr = get(handles.ms2dir,'UserData');
cd(dr{1})
cd ..
%grab it:
[fn fp] = uigetfile('*.dat','pick the time base file');
cd(currDir);
%load it:
eval(['load ' [fp fn] ' -mat'])
set(handles.timeBase,'UserData',timebaseStruc)

set(handles.instrucs,'String','time base loaded.')


% --- Executes on button press in selectSpot.
function selectSpot_Callback(hObject, eventdata, handles)
% hObject    handle to selectSpot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%tell the user whats up
set(handles.instrucs,'String','left click near spot of interest to select it.')
%get the image:
im = get(handles.slider,'UserData');
axes(handles.axes1) 
%get the image (frame) number
imagenum = get(handles.frameNum,'Value');
%get analData for spots:
%for paired spots
analData = get(handles.makeMatches,'UserData');  %this is where the matched spots struc lives. 
spotPairs = analData.spotsMat;
%spots in this frame:
spotsLogi = spotPairs(:,11) == imagenum;
%retrieve spots corresponding to this frame
currentSpots = spotPairs(spotsLogi,:); 
%here, to imitate what we do in removeAOIs, we'll make something the same
%shape and content as "aoiinfo2" & "FitData" elsewhere in the code (i.e. handles.pick
%UserData). Form: [frame average xPos yPos aoiSiaze aoiNumber]:
FitData = [currentSpots(:,11) currentSpots(:,2) currentSpots(:,5) currentSpots(:,6) currentSpots(:,2) currentSpots(:,3)]; % form [frameNum doesntMatter xPos yPos doesntMatter nucleusAssignment]
%select near the spot of interest
[a b but]=ginput(1);
% Get the aoi number for the aoi closest to where user clicked
num_closest=aoicompareFli([a b],FitData);
selectedSpotNum = num2str(num_closest);
set(handles.selectedSpotNum,'String',selectedSpotNum);
%set this number in user data and in selectedSpotNum
set(handles.selectSpot,'UserData',[num_closest imagenum]); % save bothe the nucleus number as well as the frame number
selectedSpotNum = num2str(num_closest);
set(handles.selectedSpotNum,'String',selectedSpotNum);
set(handles.instrucs,'String','Spot selected. Now enter the nucleus you''d like to assign it to and hit "Correct".');

function newNucleus_Callback(hObject, eventdata, handles)
% hObject    handle to newNucleus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of newNucleus as text
%        str2double(get(hObject,'String')) returns contents of newNucleus as a double


% --- Executes during object creation, after setting all properties.
function newNucleus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to newNucleus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in correct.
function correct_Callback(hObject, eventdata, handles)
% hObject    handle to correct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%retrieve the old and new nucleus assignments:
ogNukeV = get(handles.selectSpot,'UserData'); %this is a two member vector [nukeNum frameNum]
if ~isempty(ogNukeV)
    ogNuke = ogNukeV(1); %the old nucles assignment
    frameNum = ogNukeV(2); %frame number to adjust
else
    set(handles.instrucs,'String','No spot selected. Do that and come back.');
    return
end
newNuke = str2double(get(handles.newNucleus,'String'));
%get the analData struc
analData = get(handles.makeMatches,'UserData');
%correct the spotsMat first:
spotsMat = analData.spotsMat;
%first pick the row of both the correct frame and nucleus assignment:
logi1 = spotsMat(:,11) == frameNum & spotsMat(:,3) == ogNuke;
%reassign the nucleus for that row:
spotsMat(logi1,3) = newNuke;
%save this to analData
analData.spotsMat = spotsMat;
%Now correct the nukesMat:
nukesMat = analData.nukesMat;
%first pick the row of the correct frame for the OLD nucleus assignment:
logi2 = nukesMat(:,11) == frameNum & nukesMat(:,3) == ogNuke;
%get the spot info for this nuke at this time to reassign to the new nuke:
spotInfo = nukesMat(logi2,7:10);
%reassign the data for the OLD nucleus:
nukesMat(logi2,4) = 0; %change binary spot
nukesMat(logi2,7:10) = NaN; %make all the measures related to a spot NaN
%Next pick the row of the correct frame for the NEW nucleus assignment:
logi3 = nukesMat(:,11) == frameNum & nukesMat(:,3) == newNuke;
%reassign the data for the NEW nucleus:
nukesMat(logi3,4) = 1; %change binary spot
nukesMat(logi3,7:10) = spotInfo; %make all the measures related to the spot match that of the originally measured for the old nuke assignment
%now save this to analData
analData.nukesMat = nukesMat;
%Finally, correct cumulative interval array:
ciaStruc = makeInts(nukesMat);
analData.cia = ciaStruc.cia;

%save new analData:
set(handles.makeMatches,'UserData',analData);
%reset the image:
slider_Callback( handles.frameNum, eventdata, handles)
%update user:
set(handles.instrucs,'String','spot reassigned to a different nucleus. Continue re-assignments or save new analData.')
%return Correct spot panel values to default:
set(handles.selectedSpotNum,'String','--');
set(handles.selectSpot,'UserData',[]);

function nukeNum_Callback(hObject, eventdata, handles)
% hObject    handle to nukeNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nukeNum as text
%        str2double(get(hObject,'String')) returns contents of nukeNum as a double

% --- Executes during object creation, after setting all properties.
function nukeNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nukeNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in startTrackCorrect.
function startTrackCorrect_Callback(hObject, eventdata, handles)
% hObject    handle to startTrackCorrect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get the nuclues num input:
nukeNum = str2double(get(handles.nukeNum,'String'));
%get analData for nukes:
%for paired spots
aoifits = get(handles.trackNuclei,'UserData');  %this is where the matched spots struc lives. 
nukesData = aoifits.data;
%get the initial image (frame) number
iniImagenum = get(handles.frameNum,'Value');
%get the max frames in the time series (dont go over this)
mxFrame = get(handles.slider,'Max');
%set up a while loop for multiple frames:
posV = [];
flag = 0;
while flag == 0
    %get the image (frame) number
    imagenum = get(handles.frameNum,'Value');
    %tell the user whats up
    update = sprintf('left click where you''d like to place nucleus %d in frame %d. right click when done.',nukeNum,imagenum);
    set(handles.instrucs,'String',update)
    %get the image:
    im = get(handles.slider,'UserData');
    axes(handles.axes1) 
    %retrieve the row for teh nuke (col 1) in this frame (col 2):
    logi1 = nukesData(:,2) == imagenum & nukesData(:,1) == nukeNum;
    nukeRow = nukesData(logi1,:);
    %get the new xPos & yPos
    [xPos yPos but]=ginput(1);
    %first set up the exit, then reassign the nuke position, advance the image
    if but == 3
        confirm = sprintf('nucleus %d has been given new positions from frame %d - %d',nukeNum,iniImagenum,imagenum-1);
        set(handles.instrucs,'String',confirm);
        flag = 1;
    else
        %Reassign nuke positions (I measured the time it takes to write
        %this info within each loop or at teh end, and its the same either
        %way):
        nukeRow(4:5) = [xPos yPos];
        %insert into analData:
        try
            nukesData(logi1,:) = nukeRow;
        catch
            %avoids common errors, I believe
            set(handles.instrucs,'String','something went wrong. no re-assignments made.');
            return
        end
        aoifits.data = nukesData;
        set(handles.trackNuclei,'UserData',aoifits);
        %now advance the image:
        %advance slider value (if we're not at the end of the time series)
        if imagenum == mxFrame
            confirm = sprintf('nucleus %d has been given new positions from frame %d - %d',nukeNum,iniImagenum,imagenum);
            set(handles.instrucs,'String',confirm);
            flag = 1;
            %updated the image
            slider_Callback( handles.frameNum, eventdata, handles)
        else
            set(handles.slider,'Value',imagenum+1);
            %set displayed frame number
            display = num2str(imagenum+1);        
            set(handles.frameNum,'String',display)
            set(handles.frameNum,'Value',imagenum+1)
            %updated the image
            slider_Callback( handles.frameNum, eventdata, handles)
        end
    end
end



function aoiValue_Callback(hObject, eventdata, handles)
% hObject    handle to aoiValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of aoiValue as text
%        str2double(get(hObject,'String')) returns contents of aoiValue as a double


% --- Executes during object creation, after setting all properties.
function aoiValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to aoiValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in resetDisplay.
function resetDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to resetDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

slider_Callback( handles.frameNum, eventdata, handles)

% --- Executes on button press in addSpot.
function addSpot_Callback(hObject, eventdata, handles)
% hObject    handle to addSpot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%tell 'em how its done:
set(handles.instrucs,'String','left click on spot(s) to add spot(s), right click to save and exit')
% pixnum = str2double(get(handles.PixelNumber,'String'));
frameNum = str2double(get(handles.frameNum,'String'));
axes(handles.axes1) 
flag = 0;
%get the spots
AllSpots = get(handles.findAllSpots,'UserData');
cell = AllSpots.AllSpotsCells;
%get teh indicie pertaining to this frame:
logi = cell2mat(cell(:,3)) == frameNum;
aoinumber = cell{logi,2}+1; %number of spots in this frame (plus one)
aoiinfo = [];
while flag == 0
    [xPos yPos but]=ginput(1);
    if but==3
        flag=1;
    else
        aoiinfo=[aoiinfo; frameNum 1 xPos yPos 5 aoinumber]; %[frameNum numOfFramesAveraged xPos yPos pixelNum aoiNumber] indicies 2 & 5 dont matter
        aoinumber=aoinumber+1;                      %Give each aoi a number
        axes(handles.axes1);
        hold on 
        viscircles(aoiinfo(end,3:4),3,'Color','r','LineWidth',0.5,'EnhanceVisibility',false);
        hold off
    end
end
%get/set the fp & fn:
%fn & fp from hidden element hisDir or ms2dir:
hisImage = get(handles.hisImage,'Value');
if hisImage == 1
    fileCell = get(handles.hisDir,'UserData');
    fp = fileCell{1}; fn = fileCell{2};
    folder = [fp fn];
else
    fileCell = get(handles.ms2dir,'UserData');
    fp = fileCell{1}; fn = fileCell{2};
    folder = [fp fn];
end
handles.TiffFolder = folder;
guidata(gcbo,handles) ;
%fit the selected spots
fits = fitAspot(handles,aoiinfo); %fits format: [1.gaussAmp 2.xCenter 3.yCenter 4.gaussSigma 5.offsetFromBackground 6.integratedIntensity];
newSpotNum = size(fits,1);
%asign to amp/sigma/offset/intInt to AllSpots:
%centers:
AllSpots.AllSpotsCells{logi,1} = [cell{logi,1}; fits(:,2:3)];
%number of spots in this frame
AllSpots.AllSpotsCells{logi,2} = cell{logi,2} + newSpotNum;
%[amp sigma offset integratedInt]
AllSpots.AllSpotsCells{logi,4} = [cell{logi,4}; fits(:,[1 4 5 6])];
%update GUI handles and user
set(handles.findAllSpots,'UserData',AllSpots);
update = sprintf('spots added to frame %d',frameNum);
set(handles.instrucs,'String',update)

% --- Executes on button press in removeSpot.
function removeSpot_Callback(hObject, eventdata, handles)
% hObject    handle to removeSpot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%tell 'em how its done:
set(handles.instrucs,'String','left click near spot to remove spot, right click to save and exit')
%set up the loop:
frameNum = str2double(get(handles.frameNum,'String'));
axes(handles.axes1) 
flag = 0;
%get the spots
AllSpots = get(handles.findAllSpots,'UserData');
cell = AllSpots.AllSpotsCells;
%get teh indicie pertaining to this frame:
logi = cell2mat(cell(:,3)) == frameNum;
%number of spots in this frame
aoinumber = cell{logi,2}; 
%make something of the format FitData for this frames spots [frame ave xPos yPos size AOInum] (only position and AOInum matter}
FitData = zeros(aoinumber,6);
%add position:
FitData(:,3:4) = cell{logi,1};
%add aoi number (relative to this frame):
aoiV = [1:aoinumber];
FitData(:,6) = aoiV;
%add frame (not sure we need to do this):
FitData(:,1) = frameNum;
%run the spot remover:
while flag==0
    [a b but]=ginput(1);
    if but==3
        flag=1; %right clicking exits the loop
    else
        % Get the aoi number for the aoi closest to where user clicked
        num_closest=aoicompareFli([a b],FitData);
        % logical array, =1 when it matches the aoi number
        logik=(FitData(:,6)==num_closest); 
        % remove information for that aoi
        FitData(logik,:)=[];     
        % Update the existin list of aoi information so that no aoi numbers are skipped:
        FitData=update_FitData_aoinum(FitData);
        %for iteration decrement the aoinumber:
        aoinumber = aoinumber-1;
    end
    %update the AllSpots cell array:
    %x & y positions:
    AllSpots.AllSpotsCells{logi,1} = FitData(:,3:4);
    %number of spots in this frame is now reduced by one (in the if statement above):
    AllSpots.AllSpotsCells{logi,2} = aoinumber;
    %get rid of the gaussian fit data for spot number 'num_closest':
    if aoinumber ~=0
        %if its not the only spot in the frame:
        gData = cell{logi,4};
        gData(num_closest,:) = [];
        AllSpots.AllSpotsCells{logi,4} = gData;   
    else
        %if it is the only spot in the frame
        AllSpots.AllSpotsCells{logi,4} = [];
    end
    % Update the list of spots in the gui data:
    set(handles.findAllSpots,'UserData',AllSpots);
    % Refresh display
    slider_Callback(handles.frameNum, eventdata, handles)
end
set(handles.instrucs,'String','spots updated. Now you can save them.')


% --- Executes on button press in autoScale.
function autoScale_Callback(hObject, eventdata, handles)
% hObject    handle to autoScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autoScale


% --- Executes on button press in useTwo.
function useTwo_Callback(hObject, eventdata, handles)
% hObject    handle to useTwo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useTwo
val  = get(hObject,'Value');
% when we select this button show the useTwoPanel
if val == 1
    set(handles.useTwoPanel,'Visible','On');
elseif val == 0
    set(handles.useTwoPanel,'Visible','Off');
end

function brightnessOne_Callback(hObject, eventdata, handles)
% hObject    handle to brightnessOne (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of brightnessOne as text
%        str2double(get(hObject,'String')) returns contents of brightnessOne as a double


% --- Executes during object creation, after setting all properties.
function brightnessOne_CreateFcn(hObject, eventdata, handles)
% hObject    handle to brightnessOne (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function diameterOne_Callback(hObject, eventdata, handles)
% hObject    handle to diameterOne (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diameterOne as text
%        str2double(get(hObject,'String')) returns contents of diameterOne as a double


% --- Executes during object creation, after setting all properties.
function diameterOne_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diameterOne (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function noiseOne_Callback(hObject, eventdata, handles)
% hObject    handle to noiseOne (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noiseOne as text
%        str2double(get(hObject,'String')) returns contents of noiseOne as a double


% --- Executes during object creation, after setting all properties.
function noiseOne_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noiseOne (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function brightnessTwo_Callback(hObject, eventdata, handles)
% hObject    handle to brightnessTwo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of brightnessTwo as text
%        str2double(get(hObject,'String')) returns contents of brightnessTwo as a double


% --- Executes during object creation, after setting all properties.
function brightnessTwo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to brightnessTwo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function diameterTwo_Callback(hObject, eventdata, handles)
% hObject    handle to diameterTwo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diameterTwo as text
%        str2double(get(hObject,'String')) returns contents of diameterTwo as a double


% --- Executes during object creation, after setting all properties.
function diameterTwo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diameterTwo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function noiseTwo_Callback(hObject, eventdata, handles)
% hObject    handle to noiseTwo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noiseTwo as text
%        str2double(get(hObject,'String')) returns contents of noiseTwo as a double


% --- Executes during object creation, after setting all properties.
function noiseTwo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noiseTwo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in reverse.
function reverse_Callback(hObject, eventdata, handles)
% hObject    handle to reverse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of reverse


% --- Executes on button press in displayNoSpots.
function displayNoSpots_Callback(hObject, eventdata, handles)
% hObject    handle to displayNoSpots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of displayNoSpots


% --- Executes on button press in loadNoSpots.
function loadNoSpots_Callback(hObject, eventdata, handles)
% hObject    handle to loadNoSpots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fl = get(handles.outputText,'UserData'); % file locations are loaded at GUI launch
filestring = get(handles.input,'String');

eval(['load ' fl.data filestring ' -mat']) %this should load the struc 'AllSpots'
set(handles.loadNoSpots,'UserData',noSpotMat);

set(handles.instrucs,'String','no spots loaded');
set(handles.input,'String','default.dat')

% --- Executes on button press in intNoSpots.
function intNoSpots_Callback(hObject, eventdata, handles)
% hObject    handle to intNoSpots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% make sure there is an analData to add the output to
try 
    analData = get(handles.makeMatches,'UserData');
catch
    set(handles.instrucs,'String','no analData loaded. do that and come back.')
    return
end    
%Tell user what's up; set up a timer to display at command line
set(handles.intNoSpots,'String','...')
set(handles.instrucs,'String','integrating no spots')
pause(0.1)
%get/set the fp & fn:
%fn & fp from hidden element hisDir or ms2dir:
hisImage = get(handles.hisImage,'Value');
if hisImage == 1
    fileCell = get(handles.hisDir,'UserData');
    fp = fileCell{1}; fn = fileCell{2};
    folder = [fp fn];
else
    fileCell = get(handles.ms2dir,'UserData');
    fp = fileCell{1}; fn = fileCell{2};
    folder = [fp fn];
end
handles.TiffFolder = folder;
guidata(gcbo,handles) ;
%find the spots, now with fitting 
tic;
noSpotFitMat = AOIintegrate(handles);  
toc
set(handles.intNoSpots,'UserData',noSpotFitMat);
analData.noSpots = noSpotFitMat;
analData.noSpotsDescription = '[1.time(s) 2.frameNum 3.nucleus 4.binarySpot 5.noSpotXpos 6.noSpotYpos 7.gaussAmp 8.sigma 9.offset 10.intInt 11.relativeFrameNumber]';
set(handles.makeMatches,'UserData',analData);
%update:
set(handles.instrucs,'String','noSpots integrated. added to current analData. now resave that')
set(handles.intNoSpots,'String','Int. no spots')


% --- Executes during object creation, after setting all properties.
function pick_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
