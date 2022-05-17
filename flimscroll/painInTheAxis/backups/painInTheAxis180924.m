function varargout = painInTheAxis(varargin)
% PAININTHEAXIS MATLAB code for painInTheAxis.fig
%      PAININTHEAXIS, by itself, creates a new PAININTHEAXIS or raises the existing
%      singleton*.
%
%      H = PAININTHEAXIS returns the handle to a new PAININTHEAXIS or the handle to
%      the existing singleton*.
%
%      PAININTHEAXIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PAININTHEAXIS.M with the given input arguments.
%
%      PAININTHEAXIS('Property','Value',...) creates a new PAININTHEAXIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before painInTheAxis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to painInTheAxis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% Harden 2018
% Where shit's at:
%   directory for easy loading: handles.dir (no set/get)
%   Ant reader: get(handles.loadAntStack,'UserData')
%   Ant max projection: get(handles.maxProjAnt,'UserData'

% Edit the above text to modify the response to help painInTheAxis

% Last Modified by GUIDE v2.5 17-Sep-2018 18:11:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @painInTheAxis_OpeningFcn, ...
                   'gui_OutputFcn',  @painInTheAxis_OutputFcn, ...
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


% --- Executes just before painInTheAxis is made visible.
function painInTheAxis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to painInTheAxis (see VARARGIN)

% Choose default command line output for painInTheAxis
handles.output = hObject;

%turn off axes:
axes(handles.axes1);
axis('off');
axes(handles.axes2);
axis('off');
axes(handles.axes3);
axis('off');

%make a place to save the direcrtory where im are stored for easy
%retrieval:
handles.dir = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes painInTheAxis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = painInTheAxis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadAntStack.
function loadAntStack_Callback(hObject, eventdata, handles)
% hObject    handle to loadAntStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currDir = cd;
%if a file has already been retrieved, go to that dir
try
    dr = handles.dir;
    cd(dr{1});
end

%get the image file
[fn fp] = uigetfile('*.*','pick the ANT image');
reader = bfGetReader([fp fn]);
% im = bfGetPlane(reader, 1); %1st image

%save the dir:
handles.dir = {fp fn};
guidata(hObject, handles);

%save the reader?
set(handles.loadAntStack,'UserData',reader)

set(handles.instrucs,'String','ANT image stack loaded. Now make a max. proj.')
cd(currDir);


% --- Executes on button press in maxProjAnt.
function maxProjAnt_Callback(hObject, eventdata, handles)
% hObject    handle to maxProjAnt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

reader = get(handles.loadAntStack,'UserData');
omeMeta = reader.getMetadataStore();

% nT = omeMeta.getPixelsSizeT(0).getValue(); %number of time points
nZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
nX = omeMeta.getPixelsSizeX(0).getValue();
nY = omeMeta.getPixelsSizeY(0).getValue();
nC = omeMeta.getPixelsSizeC(0).getValue();

%for clearing heap space (needed for jheapcl.m) ... maybe it works?
% javaaddpath(which('C:\Users\tth12\matlab\fig_files\jheapcl\MatlabGarbageCollector.jar'))

%allocate:
Imax = ones(nY,nX,nC);
iT = 1; %there should only be one time point.
%max project both channels
for iC = 1:nC %one iteration takes about 4s
%     jheapcl  %right now this is throwing errors. fgure that out, then
%     uncomment both this and the one below
        fprintf('projecting channel %d \n',iC)
        %more allocation:
        I = ones(nY,nX,nZ);    
        for iZ = 1:nZ
            iPlane = reader.getIndex(iZ - 1, iC - 1, iT - 1) + 1; 
            I(:,:,iZ) = bfGetPlane(reader, iPlane);
        end
    %     end
        Imax(:,:,iC) = max(I,[],3);
    %%%%%clearing java heap space, not sure if its necessary or if it
    %%%%%works, but it doesnt hurt.
%     jheapcl
end

%plot
chanChoice = get(handles.antChan1,'Value');
if chanChoice == 1
    chan = 1;
else
    chan = 2;
end
im = Imax(:,:,chan);
hAx = handles.axes1;
imshow(imadjust(mat2gray(im,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx); %this works and avoid having to fuck with scale bar thingies 
%alternatively:
%     axes(handles.axes1);
%     imagesc(im,[0 255] );axis('equal');axis('off'); colormap(gray(256));

%up to here works well. figure out saving. add channel option. then loading
%and showing. pretty easy so far.

%save the max proj image in the GUI:
set(handles.maxProjAnt,'UserData',Imax)

%save the image to disk:
fp = handles.dir{1};
fn = extractBefore(handles.dir{2},'_'); %try this, should work well enough to get some info on the image to the filename    
fileOutPath = [fp fn '-antMaxProj.tiff'];
bfsave(Imax,fileOutPath, 'BigTiff', true);

set(handles.instrucs,'String','ANT max proj created and saved to disk in the same dir as the image stack')


% --- Executes on button press in plotAnt.
function plotAnt_Callback(hObject, eventdata, handles)
% hObject    handle to plotAnt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Imax = get(handles.maxProjAnt,'UserData');
if isempty(Imax)
    set(handles.instrucs,'String','no max proj image. Make or load one')
    return
end

chanChoice = get(handles.antChan1,'Value');
if chanChoice == 1
    chan = 1;
else
    chan = 2;
end
im = Imax(:,:,chan);
hAx = handles.axes1;
imshow(imadjust(mat2gray(im,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);

% --- Executes on button press in loadAntMax.
function loadAntMax_Callback(hObject, eventdata, handles)
% hObject    handle to loadAntMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currDir = cd;
%if a file has already been retrieved, go to that dir
try
    dr = handles.dir;
    cd(dr{1});
end

%get the image file
[fn fp] = uigetfile('*.*','pick the ANT max projection');
% eval(['load ' fp fn])
% eval(['load ' foldstruc.gfolder 'header.mat'])
reader = bfGetReader([fp fn]);
%contained in the reader is two images corresponding to the max proj in
%each channel. grab them both:
im1 = bfGetPlane(reader, 1); %1st image
im2 = bfGetPlane(reader, 2);

%now put them in the same mat:
Imax(:,:,1) = im1;
Imax(:,:,2) = im2;

%save the dir:
handles.dir = {fp fn};
guidata(hObject, handles);

%save the image to GUI data
set(handles.maxProjAnt,'UserData',Imax)

%show the image:
chanChoice = get(handles.antChan1,'Value');
if chanChoice == 1
    chan = 1;
else
    chan = 2;
end
im = Imax(:,:,chan);
hAx = handles.axes1;
imshow(imadjust(mat2gray(im,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);

set(handles.instrucs,'String','ANT max projection loaded.')
cd(currDir);

% --- Executes on button press in loadPostStack.
function loadPostStack_Callback(hObject, eventdata, handles)
% hObject    handle to loadPostStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currDir = cd;
%if a file has already been retrieved, go to that dir
try
    dr = handles.dir;
    cd(dr{1});
end

%get the image file
[fn fp] = uigetfile('*.*','pick the POST image');
reader = bfGetReader([fp fn]);
% im = bfGetPlane(reader, 1); %1st image

%save the dir:
handles.dir = {fp fn};
guidata(hObject, handles);

%save the reader?
set(handles.loadPostStack,'UserData',reader)

set(handles.instrucs,'String','POST image stack loaded. Now make a max. proj.')
cd(currDir);

% --- Executes on button press in maxProjPost.
function maxProjPost_Callback(hObject, eventdata, handles)
% hObject    handle to maxProjPost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

reader = get(handles.loadPostStack,'UserData');
omeMeta = reader.getMetadataStore();

% nT = omeMeta.getPixelsSizeT(0).getValue(); %number of time points
nZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
nX = omeMeta.getPixelsSizeX(0).getValue();
nY = omeMeta.getPixelsSizeY(0).getValue();
nC = omeMeta.getPixelsSizeC(0).getValue();

%for clearing heap space (needed for jheapcl.m) ... maybe it works?
% javaaddpath(which('C:\Users\tth12\matlab\fig_files\jheapcl\MatlabGarbageCollector.jar'))

%allocate:
Imax = ones(nY,nX,nC);
iT = 1; %there should only be one time point.
%max project both channels
for iC = 1:nC %one iteration takes about 4s
%     jheapcl  %right now this is throwing errors. fgure that out, then
%     uncomment both this and the one below
        fprintf('projecting POST channel %d \n',iC)
        %more allocation:
        I = ones(nY,nX,nZ);    
        for iZ = 1:nZ
            iPlane = reader.getIndex(iZ - 1, iC - 1, iT - 1) + 1; 
            I(:,:,iZ) = bfGetPlane(reader, iPlane);
        end
    %     end
        Imax(:,:,iC) = max(I,[],3);
    %%%%%clearing java heap space, not sure if its necessary or if it
    %%%%%works, but it doesnt hurt.
%     jheapcl
end

%plot
chanChoice = get(handles.postChan1,'Value');
if chanChoice == 1
    chan = 1;
else
    chan = 2;
end
im = Imax(:,:,chan);
hAx = handles.axes2;
imshow(imadjust(mat2gray(im,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx); %this works and avoid having to fuck with scale bar thingies 
%alternatively:
%     axes(handles.axes1);
%     imagesc(im,[0 255] );axis('equal');axis('off'); colormap(gray(256));

%up to here works well. figure out saving. add channel option. then loading
%and showing. pretty easy so far.

%save the max proj image in the GUI:
set(handles.maxProjPost,'UserData',Imax)

%save the image to disk:
fp = handles.dir{1};
fn = extractBefore(handles.dir{2},'_'); %try this, should work well enough to get some info on the image to the filename    
fileOutPath = [fp fn '-postMaxProj.tiff'];
bfsave(Imax,fileOutPath, 'BigTiff', true);

set(handles.instrucs,'String','POST max proj created and saved to disk in the same dir as the image stack')

% --- Executes on button press in plotPost.
function plotPost_Callback(hObject, eventdata, handles)
% hObject    handle to plotPost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Imax = get(handles.maxProjPost,'UserData');
if isempty(Imax)
    set(handles.instrucs,'String','no max proj image. Make or load one')
    return
end

chanChoice = get(handles.postChan1,'Value');
if chanChoice == 1
    chan = 1;
else
    chan = 2;
end
im = Imax(:,:,chan);
hAx = handles.axes2;
imshow(imadjust(mat2gray(im,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);

% --- Executes on button press in loadPostMax.
function loadPostMax_Callback(hObject, eventdata, handles)
% hObject    handle to loadPostMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currDir = cd;
%if a file has already been retrieved, go to that dir
try
    dr = handles.dir;
    cd(dr{1});
end

%get the image file
[fn fp] = uigetfile('*.*','pick the POST max projection');
% eval(['load ' fp fn])
% eval(['load ' foldstruc.gfolder 'header.mat'])
reader = bfGetReader([fp fn]);
%contained in the reader is two images corresponding to the max proj in
%each channel. grab them both:
im1 = bfGetPlane(reader, 1); %1st image
im2 = bfGetPlane(reader, 2);

%now put them in the same mat:
Imax(:,:,1) = im1;
Imax(:,:,2) = im2;

%save the dir:
handles.dir = {fp fn};
guidata(hObject, handles);

%save the image to GUI data
set(handles.maxProjPost,'UserData',Imax)

%show the image:
chanChoice = get(handles.postChan1,'Value');
if chanChoice == 1
    chan = 1;
else
    chan = 2;
end
im = Imax(:,:,chan);
hAx = handles.axes2;
imshow(imadjust(mat2gray(im,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);

set(handles.instrucs,'String','POST max projection loaded.')
cd(currDir);
