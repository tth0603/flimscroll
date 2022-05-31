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
%
%   Bodies buried:
%       1. cannot track through mitosis. each n.c. must be analyzed
%       individually. --will likely never be addressed
%       2. alignment through the painInTheAxis GUI is still a bit in beta.
%       3. Only can handle max proj data. will need to be improved. 
%       4. Right now we manually input a time base (ttb). perhaps we can
%       improve? --fixed
%       5. N.C. is also not considered in the the ttb - add NC frame inputs
%       and NC to analData --fixed
%       6. need a way to evaluate the DL. add this to nuclei box --fixed
%       7. get rid of nuclei doublets
%       
% Timothy Harden 2018
% 

% Last Modified by GUIDE v2.5 10-Oct-2018 14:57:16

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
load('fileLocations.dat','-mat','fileLocations')
set(handles.outputText,'UserData',fileLocations);

%turn off axes:
axes(handles.axes1);
axis('off');

%initialize some stuff that is not is handles. esp good for fitting
% handles.FitData = [];
handles.RollingBallRadius = 15;
handles.RollingBallHeight = 5; %if Choice is a thing, then we can also switch these
handles.Pixnums = [];
% handles.TrackAOIs = 1; %this is now the tag for the hidden text titled "Fit Choice"
handles.TiffFolder = '';

version = 1;
fprintf('version %d \n',version);

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
aoi = get(handles.aoiValue,'Value'); %now [xmin xmax ymin ymax]

%decide if it's MS2 or his:
hisVal = get(handles.hisImage,'Value');
ms2Val = get(handles.ms2image,'Value');

%get the image scale stuff
clowval = get(handles.minSlider,'Value');
chival = get(handles.maxSlider,'Value');

%show numbers?
shoNum = get(handles.aoiNum,'Value');

%feature display size
pixnum = str2double(get(handles.PixelNumber,'String')); %how big of spots

if ms2Val
    %get bf reader:
    reader = get(handles.load,'UserData');
    %get image

%    im = bfGetPlane(reader,reader.getIndex(0, 1, imagenum - 1));  %in the past this has worked. I am growing a bit tired of how finicky this bf shit is
    im = bfGetPlane(reader,reader.getIndex(imagenum - 1, 0, 0)+1);       %getIndex is (Z, C, T) I think. the channel should be 0 for gfp, but may change
    %plot
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

elseif hisVal
        %get bf reader:
        reader = get(handles.load2,'UserData');
        %get image
%         im = bfGetPlane(reader,reader.getIndex(0, 1, imagenum - 1));   
        im = bfGetPlane(reader,reader.getIndex(imagenum - 1, 0, 0)+1);
        %plot
        hAx = handles.axes1;
        if crp && length(aoi) == 4
            %smlIm = im(aoi(2):aoi(4)+aoi(2),aoi(1):aoi(3)+aoi(1)); %weird format of aoi bc the axis of the image is inverted, i think
            %imshow(imadjust(mat2gray(smlIm,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);
            imagesc(im,[clowval chival] );axis('equal');axis('off');colormap(gray(256));axis(aoi);  %axis([xmin xmax ymin ymax])
        else
            %imshow(imadjust(mat2gray(im,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);
            imagesc(im,[clowval chival] );axis('equal');axis('off');%colormap(gray(256));%axis(limitsxy)
        end
        set(handles.slider,'UserData',im);

else
    set(handles.instrucs,'String','something went wrong with MS2/his toggle')
end

% %if we have spots from feature picker, plot themL
% %heres where feature spots live:
% featureSpots = get(handles.pick,'UserData');
% %plotting:
% if ~isempty(featureSpots)
%     aoirose = size(featureSpots,1);
%     for indx=1:aoirose
%         viscircles(featureSpots(indx,3:4),pixnum/2,'Color','b','LineWidth',0.5,'EnhanceVisibility',false);
%         if shoNum == 1
%             %showing numbers not yet implemented here
%         end
%     end 
% end

%show tracked nuclei:
if get(handles.displayTracks,'Value') == 1
    aoifits = get(handles.trackNuclei,'UserData');
    if ~isempty(aoifits)
        aoirose = length(aoifits.aoiinfo2);
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
        %this is where the matched spots struc lives. current form: [1.time 2.frameNum 3.nucleus 4.distanceFromNucleus 5.xPos 6.yPos 7.amplitude 8.sigma 9.offset 10.intensity 11.relativeFrameNumber 12.APpos]
        spotPairs = analData.spotsMat;
        spotsLogi = spotPairs(:,2) == imagenum;
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

aoi = get(handles.aoiValue,'Value'); %[XMIN YMIN WIDTH HEIGHT]


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
newVal = num2str(val - 10);
set(handles.brightness,'String',newVal);
%re-pick:
pick_Callback(hObject, eventdata,handles);

% --- Executes on button press in brightnessPlus10.
function brightnessPlus10_Callback(hObject, eventdata, handles)
% hObject    handle to brightnessPlus10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = str2double(get(handles.brightness,'String'));
newVal = num2str(val + 10);
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
im = get(handles.slider,'UserData');

axes(handles.axes1) 
flag = 0;
aoiinfo = [];

%get the aois, if none picked, set instrucs:
FitData = get(handles.pick,'UserData'); %this may be an empty struc. 

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

pixnum = str2double(get(handles.PixelNumber,'String'));
frameNum = str2double(get(handles.frameNum,'String'));
ave = 1; %use this later for frame averageing
%get the image:
im = get(handles.slider,'UserData');

axes(handles.axes1) 
flag = 0;
aoiinfo = [];

%get the aois, if none picked, set instrucs:
FitData = get(handles.pick,'UserData'); %this may be an empty struc. 

while flag==0
    [a b but]=ginput(1);
    if but==3
        flag=1;
    else
                                            % Get the aoi number for the
                                            % aoi closest to where user
                                            % clicked
       
        num_closest=aoicompareFli([a b],FitData);
       
       
                                            % logical array, =1 when it
                                            % matches the aoi number
       
        logik=(FitData(:,6)==num_closest);
       
        FitData(logik,:)=[];          % remove information for that aoi
                                            % Refresh display
      
        FitData=update_FitData_aoinum(FitData);
        
        %slider1_Callback(handles.ImageNumber, eventdata, handles, dum,images,folder)
%         slider_Callback(handles.frameNum, eventdata, handles)
        
    end
FitData=[FitData; aoiinfo];        %TH not sure why we do this
                                                    % Update the existin list of aoi 
                                                  % information so that no
                                                  % aoi numbers are skipped

FitData=update_FitData_aoinum(FitData);
set(handles.pick,'UserData',FitData);
slider_Callback(handles.frameNum, eventdata, handles)

set(handles.instrucs,'String','AOI list updated')

end


% --- Executes on button press in findAllSpots.
function findAllSpots_Callback(hObject, eventdata, handles)
% hObject    handle to findAllSpots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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
AllSpots=FindAllSpotsFli(handles,3500);  
set(handles.findAllSpots,'UserData',AllSpots);
toc

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

eval(['load ' fl.data filestring ' -mat']) %this should load the struc 'AllSpots'
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
set(handles.instrucs,'String','nuclei locations saved');
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
set(handles.instrucs,'String','tracking. this may take a moment. 0.0014s /frame /pix /nuclei on my machine')
pause(0.1)

%some inputs:
pixnum = str2double(get(handles.PixelNumber,'String'));
frms = eval([get(handles.FrameRange,'String') ';']);
ave = 1; %use this later for frame averageing
aoiinf = get(handles.pick,'UserData');
FitData = get(handles.pick,'UserData');

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
[mfrms nfrms]=size(frms); 
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

% Build a 2D mapstruc to direct data processing, also check for existance
% of DL:
useDL = get(handles.useDriftlist,'Value');
if useDL == 0
    mapstruc2d=build_2d_mapstruc_aois_frmsFli(folder,pixnum,frms,ave,FitData,1,handles);      %no DL  
else
    %make sure we have a DL
    DL = get(handles.createDriftlist,'UserData');
    if isempty(DL)
        set(handles.instrucs,'String','No drift list loaded. Make or load that then try again')
        return;
    else
        mapstruc2d=build_2d_mapstruc_aois_frmsFli(folder,pixnum,frms,ave,FitData,2,handles);        % use a measure of global nuclei drift (a DL) to track aois. its mo betta
        %all the DL magic in the tracking pipeline takes place in the
        %sub-script build_mapstruc_cell_columnFli
    end 
end

DataOutput2d=gauss2d_mapstruc2d_v2Fli(mapstruc2d,handles); % Process the data (integrate, fit etc)
                                                   % V.2 is parallel processing  

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

eval(['load ' fl.data filestring ' -mat'])

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


% --- Executes on button press in DecrementRollingBallRadius.
function DecrementRollingBallRadius_Callback(hObject, eventdata, handles)
% hObject    handle to DecrementRollingBallRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in IncrementRollingBallHeight.
function IncrementRollingBallHeight_Callback(hObject, eventdata, handles)
% hObject    handle to IncrementRollingBallHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in DecrementRollingBallHeight.
function DecrementRollingBallHeight_Callback(hObject, eventdata, handles)
% hObject    handle to DecrementRollingBallHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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

eval(['load ' fl.data filestring ' -mat'])

%what are we calling the match output? ...store it in makeMatches user data
set(handles.makeMatches,'UserData',analData);

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
%plot new nuclei
slider_Callback( handles.frameNum, eventdata, handles) 

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
vid.ttb = get(handles.timeBase,'UserData'); %originally ttb was a row vector in some absolute time (ms) so you have to subtract the value of the first mamber. in flimscroll, its also a col vector, in units of (s), with whatever the user puts as time zero. 
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
DL = drifts.diffdriftlist;

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
spotsMat(:,12)=APpositions./APlength; %check for sign errors in the future. they will likely arise from 'angles'. soln uis likely just to make this 'abs(... '
%repeat for nuclei:
nukeAngles=atan2((-nukesMat(:,6)+coordA(2)),(nukesMat(:,5)-coordA(1)));
%Distance between the points and the ANT point:
nukeDistances=sqrt((coordA(2)-nukesMat(:,6)).^2+(coordA(1)-nukesMat(:,5)).^2);
%get dist along AP axis:
nukeAPpositions=nukeDistances.*cos(nukeAngles-APangle);
%add to mat:
nukesMat(:,12)=nukeAPpositions./APlength; %check for sign errors in the future. they will likely arise from 'angles'. soln uis likely just to make this 'abs(... '

%update descriptions:
analData.spotsMatDescription = '[1.time 2.frameNum 3.nucleus 4.distanceFromNucleusCntr 5.xPos 6.yPos 7.gaussAmp 8.sigma 9.offset 10.integratedIntensity 11.relativeFrameNumber 12.APpos]';
analData.nukesMatDescription = '[1.time 2.frameNum 3.nucleus 4.binarySpot 5.nukeXpos 6.nukeYpos 7.spotGaussAmp 8.spotSigma 9.spotOffset 10.spotIntegratedIntensity 11.relativeFrameNumber 12.nuclearAPpos]';
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
fn = extractBefore(handles.dir{2},'_'); %try this, should work well enough to get some info on the image to the filename    
fileOutPath = [fp fn '-ttb.dat'];
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
eval('load ' fp fn ' -mat')
set(handles.timeBase,'UserData',timebaseStruc)

set(handles.instruc,'String','time base loaded.')








