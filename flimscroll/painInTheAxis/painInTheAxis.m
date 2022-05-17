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
%   im stack reader: get(handles.loadStack,'UserData')
%   Ant max projection: handles.antProjection (not in traditional GUI data)
%   Post max projection: handles.postProjection (not in traditional GUI data)
%   last t point max proj: handles.lastT (not in trad GUI data)
%   hand selected points for alignment: handles.antPoint; handles.postPoint; handles.lastTpoint (not in trad GUI data)
%   AP alignment info: handles.align,'UserData' (a struc named align with two images and two positional thingys)
%   A & P position (in ANT image coords): handles.pickAP,'UserData' (a struc with two members: A and P position vectors.)
%   A & P position (in data image coords, i.e. the output): handles.save,'UserData'

% Edit the above text to modify the response to help painInTheAxis

% Last Modified by GUIDE v2.5 17-Jul-2019 16:18:08

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
% axes(handles.axes1);
% axis('off');
% axes(handles.axes2);
% axis('off');

%make a place to save the direcrtory where im are stored for easy
%retrieval:
handles.dir = [];

%make places for images & stuff in GUI handles:
handles.antProjection = [];
handles.postProjection = [];
handles.lastT = [];
handles.antPoint = []; 
handles.postPoint = [];
handles.lastTpoint = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes painInTheAxis wait for user response (see UIRESUME)
uiwait(handles.figure1);  %uncommented this to pass output to flimscroll

% --- Outputs from this function are returned to the command line.
function varargout = painInTheAxis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%to pass coord struc to flimscroll:
%default
varargout{1} = handles.output;
varargout{1} = get(handles.save,'UserData');
% delete(handles.figure1);
delete(hObject);

%got this from https://blogs.mathworks.com/videos/2010/02/12/advanced-getting-an-output-from-a-guide-gui/

% --- Executes on button press in loadStack.
function loadStack_Callback(hObject, eventdata, handles)
% hObject    handle to loadStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currDir = cd;
%if a file has already been retrieved, go to that dir
try
    dr = handles.dir;
    cd(dr{1});
    cd ..  %go up one dir for generality?
end

switch get(handles.imStackMenu,'Value')
    case 1
        %get the image file
        [fn fp] = uigetfile('*.*','pick the ANT image');
        reader = bfGetReader([fp fn]);

        %save the dir:
        handles.dir = {fp fn};
        guidata(hObject, handles);

        %save the reader?
        set(handles.loadStack,'UserData',reader)

        set(handles.instrucs,'String','ANT image stack loaded. Now make a max. proj.')
    case 2
        %get the image file
        [fn fp] = uigetfile('*.*','pick the POST image');
        reader = bfGetReader([fp fn]);
        % im = bfGetPlane(reader, 1); %1st image

        %save the dir:
        handles.dir = {fp fn};
        guidata(hObject, handles);

        %save the reader?
        set(handles.loadStack,'UserData',reader)

        set(handles.instrucs,'String','POST image stack loaded. Now make a max. proj.')
end
%make the menu disappear so the user doesn;t change it. this is because
%both the ant and post image stack readers are saved to the same place in
%the GUI data, so we don't want to make things confusing. 
% set(handles.imStackMenu,'Visible','Off') %changed my mind 190716
set(handles.maxProj,'Visible','On')

cd(currDir);

% --- Executes on button press in maxProj.
function maxProj_Callback(hObject, eventdata, handles)
% hObject    handle to maxProj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%update the user
set(handles.instrucs,'String','Making maximum projection');
pause(0.1)
%grab the meta data
reader = get(handles.loadStack,'UserData');
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
chanChoice = get(handles.chan1,'Value');
if chanChoice == 1
    chan = 1;
else
    chan = 2;
end
im = Imax(:,:,chan);
hAx = handles.axes1;
clowval = min(min(im));
chival = max(max(im));
axes(hAx);
imagesc(im,[clowval chival]);axis(hAx,'equal','off');colormap(gray(256));
% imshow(imadjust(mat2gray(im,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx); %this works and avoid having to fuck with scale bar thingies 
%alternatively:
%     axes(handles.axes1);
%     imagesc(im,[0 255] );axis('equal');axis('off'); colormap(gray(256));

switch get(handles.imStackMenu,'Value')
    case 1
        %save the max proj image in the GUI:
        %set(handles.maxProj,'UserData',Imax)
        %get rid of gross outliers for im scaling and display purposes.
        %This has a likelihood of causing headaches in the future!
        Imax(Imax > 1000) = 0;
        %save to gui data
        handles.antProjection = Imax;
        %save the image to disk:
        fp = handles.dir{1};
        fn = extractBefore(handles.dir{2},'_'); %try this, should work well enough to get some info on the image to the filename    
        fileOutPath = [fp fn '-antMaxProj.tiff'];
        bfsave(Imax,fileOutPath, 'BigTiff', true);

        set(handles.instrucs,'String','ANT max proj created and saved to disk in the same dir as the image stack')
    case 2
        %get rid of gross outliers for display purposes. Has high
        %likelihood of bugging out in the future! THis is at three other
        %instances in this gui
        Imax(Imax > 1000) = 0;
        %save the max proj image in the GUI:
        handles.postProjection = Imax;
        %save the image to disk:
        fp = handles.dir{1};
        fn = extractBefore(handles.dir{2},'_'); %try this, should work well enough to get some info on the image to the filename    
        fileOutPath = [fp fn '-postMaxProj.tiff'];
        bfsave(Imax,fileOutPath, 'BigTiff', true);

        set(handles.instrucs,'String','POST max proj created and saved to disk in the same dir as the image stack')
end

guidata(hObject, handles);
set(handles.maxProj,'Visible','Off')
set(handles.imStackMenu,'Visible','On')



% --- Executes on button press in display.
function display_Callback(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hold off

%where to display:
% if get(handles.display2,'Value')
%     hAx = handles.axes2;
% else 
%     hAx = handles.axes1;
% end
%190716 just display in 1. made display 2 radio button a vestige:
hAx = handles.axes1;

switch get(handles.displayMenu,'Value')
    case 1 %anterior
        Imax = handles.antProjection;
        if isempty(Imax)
            set(handles.instrucs,'String','no max proj image. Make or load one')
            return
        end
        %show selected channel:
        chanChoice = get(handles.chan1,'Value');
        if chanChoice == 1
            chan = 1;
        else
            chan = 2;
        end
        im = Imax(:,:,chan);
        %hAx = handles.axes1;
        clowval = min(min(im));
        chival = max(max(im));
        axes(hAx);
        imagesc(im,[clowval chival] );axis('equal');axis('off');colormap(gray(256));
%         imshow(imadjust(mat2gray(im),[0 64000]),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);

        %this is some BS attempt to get better alignment: (see
        %https://www.mathworks.com/matlabcentral/answers/89960-image-alignment-code-in-matlab)
%         stats = regionprops('table',im,'Centroid', ...
%                          'MajorAxisLength','MinorAxisLength');
% 
%       % Get centers and radii of the circles
%       centers = stats.Centroid;
%       diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
%       radii = diameters/2;
% 
%       % Plot the circles
%       hold on
%       viscircles(centers,radii);
%       hold off      
      
    case 2 %posterior
        Imax = handles.postProjection;
        if isempty(Imax)
            set(handles.instrucs,'String','no max proj image. Make or load one')
            return
        end

        chanChoice = get(handles.chan1,'Value');
        if chanChoice == 1
            chan = 1;
        else
            chan = 2;
        end
        im = Imax(:,:,chan);
        clowval = min(min(im));
        chival = max(max(im));
        axes(hAx);
        imagesc(im,[clowval chival] );axis('equal');axis('off');colormap(gray(256));
%         imshow(imadjust(mat2gray(im,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);
        
        %some bullshit, see case 2 above:
%                 stats = regionprops('table',im,'WeightedCentroid', ...
%                          'MajorAxisLength','MinorAxisLength');
% 
%       % Get centers and radii of the circles
%       centers = stats.Centroid;
%       diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
%       radii = diameters/2;
% 
%       % Plot the circles
%       hold on
%       viscircles(centers,radii,'b');
%       hold off  
      
    case 3 %last T
        Imax = handles.lastT;
        if isempty(Imax)
            set(handles.instrucs,'String','no max proj image. Make or load one')
            return
        end
        %hAx = handles.axes1;
        %im = Imax;
        clowval = min(min(Imax));
        chival = max(max(Imax));
        axes(hAx);
        imagesc(Imax,[clowval chival] );axis('equal');axis('off');colormap(gray(256));
%         imshow(imadjust(mat2gray(Imax,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);
        
    case 4 %A-P alignment
        %get the aligned images:
        align  = get(handles.align,'UserData');
        if isempty(align)
            set(handles.instrucs,'String','no aligment info. Make one first.')
            return
        end
        pTranslated = align.pTranslated;
        pTranslatedRef = align.pTranslatedRef;
        aIm = align.aIm;
        RA = align.RA;
        %show it:
        %hAx = handles.axes1;
        imshowpair(pTranslated, pTranslatedRef, aIm, RA,'Scaling','joint','Parent',hAx);shg
%         axis('off');
        %add AP axis if it exists:
        AP = get(handles.pickAP,'UserData');
        if ~isempty(AP)
            hold on;
            xs = [AP.coordA(1) AP.coordP(1)];
            ys = [AP.coordA(2) AP.coordP(2)];
            %labeling the selections:
            txt1 = ['ANT'];
            txt2 = ['POST'];
            text(xs(1),ys(1),txt1);
            text(xs(2),ys(2),txt2);
            plot(xs(1:2), ys(1:2), '-+r','LineWidth',1,'MarkerSize',12);shg %,'Parent',hAx
        end
        
        set(handles.instrucs,'String','AP alignment displayed')
        
    case 5 %Ant-last t point alignment
        %get the aligned images:
        align  = get(handles.align,'UserData');        
        if isempty(align)
            set(handles.instrucs,'String','no aligment info. Make one first.')
            return
        end
        datTranslated = align.datTranslated;
        datTranslatedRef = align.datTranslatedRef;
        aIm = align.aIm;
        RA = align.RA;
        %show it:
        %hAx = handles.axes1;
        imshowpair(datTranslated, datTranslatedRef, aIm, RA,'Scaling','joint','Parent',hAx);shg
%         axis('off');
        set(handles.instrucs,'String','AP-data alignment displayed')        
        
end

% if get(handles.display2,'Value')
%     hAx = handles.axes2;
%     imshow(imadjust(mat2gray(im,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);
% end

% --- Executes on button press in loadMax.
function loadMax_Callback(hObject, eventdata, handles)
% hObject    handle to loadMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currDir = cd;
%if a file has already been retrieved, go to that dir
try
    dr = handles.dir;
    cd(dr{1});
    cd .. %go to one dir above where the ant(or post) file lives
end
hold off
switch get(handles.displayMenu,'Value')
    case 1 %ant
        %get the image file
        [fn fp] = uigetfile('*.*','pick the ANT max projection');
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
        %save the image to GUI data
        %set(handles.maxProj,'UserData',Imax)
        %get rid of gross outliers for im scaling and display purposes.
        %This has a likelihood of causing headaches in the future!
        Imax(Imax > 1000) = 0;
        %save to gui data
        handles.antProjection = Imax;
        %show the image:
        chanChoice = get(handles.chan1,'Value');
        if chanChoice == 1
            chan = 1;
        else
            chan = 2;
        end
        im = Imax(:,:,chan);
        hAx = handles.axes1;
        clowval = min(min(im));
        chival = max(max(im));
        axes(hAx);
        imagesc(im,[clowval chival]);axis(hAx,'equal','off');colormap(gray(256));
%         imshow(imadjust(mat2gray(im,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);

        set(handles.instrucs,'String','ANT max projection loaded.')
    case 2 %post
        %get the image file
        [fn fp] = uigetfile('*.*','pick the POST max projection');
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
        %get rid of gross outliers for display purposes. Has high
        %likelihood of bugging out in the future! THis is at three other
        %instances in this gui
        Imax(Imax > 1000) = 0;
        %save the image to GUI data
        %set(handles.maxProjPost,'UserData',Imax)
        handles.postProjection = Imax;
        %show the image:
        chanChoice = get(handles.chan1,'Value');
        if chanChoice == 1
            chan = 1;
        else
            chan = 2;
        end
        im = Imax(:,:,chan);
        hAx = handles.axes1;
        clowval = min(min(im));
        chival = max(max(im));
        axes(hAx);
        imagesc(im,[clowval chival]);axis(hAx,'equal','off');colormap(gray(256));
%         imshow(imadjust(mat2gray(im,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);

        set(handles.instrucs,'String','POST max projection loaded.')
    case 3 %last T point
        %this one is a bit different. this max proj already exists in
        %previously projected image series, so pick that file and go grab
        %the last image
        set(handles.instrucs,'String','pick the max proj image series youd like to use')
        pause(0.1)
        
        [fn fp] = uigetfile('*.*','pick the max proj image series youd like to use');
        reader = bfGetReader([fp fn]);
        meta4 = reader.getMetadataStore();
        nT = meta4.getPixelsSizeZ(0).getValue(); %for some reason maxProjLS is written such that the t points are labeled z stacks, so here we adjust accordingly. I am sure I had a good reason for doing it this way?
        iMax = bfGetPlane(reader, nT);
        %save it:
        handles.lastT = iMax;
        %show it
        axes(handles.axes1);
        hAx = handles.axes1;
        clowval = min(min(iMax));
        chival = max(max(iMax));
        axes(hAx);
        imagesc(iMax,[clowval chival]);axis(hAx,'equal','off');colormap(gray(256));
%         imshow(imadjust(mat2gray(iMax,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);
        set(handles.instrucs,'String','Last T max proj retrieved.')
        
    case 4 %AP alignment
        set(handles.instrucs,'String','this fn is not yet supported. bug harden')
        
    case 5 %AP-data alignment
        set(handles.instrucs,'String','this fn is not yet supported. bug harden')
        
end
guidata(hObject, handles);
cd(currDir);


% --- Executes on selection change in displayMenu.
function displayMenu_Callback(hObject, eventdata, handles)
% hObject    handle to displayMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns displayMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from displayMenu


% --- Executes during object creation, after setting all properties.
function displayMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to displayMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in imStackMenu.
function imStackMenu_Callback(hObject, eventdata, handles)
% hObject    handle to imStackMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns imStackMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imStackMenu


% --- Executes during object creation, after setting all properties.
function imStackMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imStackMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in align.
function align_Callback(hObject, eventdata, handles)
% hObject    handle to align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%seed alignment?
% no seeding, manual alignment:
if get(handles.seed,'Value') == 0 && get(handles.manual,'Value') == 0    

    switch get(handles.alignMenu,'Value')
        case 1 %A-P
            tic;
            set(handles.instrucs,'String','aligning ANT & POST images. this may take a minute')
            set(handles.align,'String','. . .')
            pause(0.1)

            %what channel:
            chanChoice = get(handles.chan1,'Value');
            if chanChoice == 1
                chan = 1;
            else 
                chan = 2;
            end
            %images:
            aIm = handles.antProjection(:,:,chan);
            pIm = handles.postProjection(:,:,chan);
            %rescale images for display:
            pIm = imadjust(mat2gray(pIm,[0 64000]));
            aIm = imadjust(mat2gray(aIm,[0 64000]));
            %from findAPaxisCSv11 that I wrote previously:
            %TH set up alignment params:
            [optimizer, metric] = imregconfig('multimodal');
            %here we reduce the InitialRadius property of the optimizer to ensure we
            %get an answer that converges
            r = optimizer.InitialRadius;
            optimizer.InitialRadius = r*10^(-1);
            optimizer.MaximumIterations = 500;
    %         optimizer.GrowthFactor = 1.5;
    %         optimizer.Epsilon = 1.5*10^(-10);
            %try messing with metric:
    %         metric.NumberOfSpatialSamples = 50000;
    %         metric.NumberOfHistogramBins = 500;
            % this gives you the transformation matrix
            tform = imregtform(pIm, aIm, 'translation', optimizer, metric);       %syntax: imregtform(movingPic, referencePic, ... )
            %shift the POST image:
            xShift = tform.T(3,1);
            yShift = tform.T(3,2);      
            [pTranslated,pTranslatedRef] = imtranslate(pIm,[xShift yShift],'OutputView','full');
            %TH pTranslated is the translated image; pTranslatedRef is an imref2D
            %object to tell MatLab where to place the image in a figure.
            RA = imref2d(size(aIm)); %TH creates the inref2D object for aMaxProj image. each of these objects are used by imshowpair, amung other fns, to place images correctly in figures
            %Plot the alignment, select the ANT and POST position in ref frame of zoomed out, anterior image:
            hAx = handles.axes1;
%             axis('off');
            imshowpair(pTranslated, pTranslatedRef,aIm,RA,'Scaling','independent','Parent',hAx);shg
            %store:
            align.pTranslated = pTranslated;
            align.pTranslatedRef = pTranslatedRef;
            align.aIm = aIm;
            align.RA = RA;            
            set(handles.align,'UserData',align) 
            toc

            set(handles.instrucs,'String','images aligned. now pick the ANT & POST positions. Did alignment not work? try seeding it by selecting a corresponding point in Ant & Post images.')
            set(handles.align,'String','Align')

        case 2 %Ant last T point
            tic;
            set(handles.instrucs,'String','aligning the last T point image to A-P image')
            set(handles.align,'String','. . .')
            pause(0.1)
            %what channel:
            chanChoice = get(handles.chan1,'Value');
            if chanChoice == 1
                chan = 1;
            else 
                chan = 2;
            end
            aIm = handles.antProjection(:,:,chan);
            datIm = handles.lastT; %this image onlyu has one channel
            %rescale images for display:
            datIm = imadjust(mat2gray(datIm,[0 64000]));
            aIm = imadjust(mat2gray(aIm,[0 64000]));
            %TH set up alignment params:
            [optimizer, metric] = imregconfig('multimodal');
            %here we reduce the InitialRadius property of the optimizer to ensure we
            %get an answer that converges
            r = optimizer.InitialRadius;
            optimizer.InitialRadius = r*10^(-1);
            optimizer.MaximumIterations = 500;
            % this gives you the transformation matrix
            tform = imregtform(datIm, aIm, 'translation', optimizer, metric);       %syntax: imregtform(movingPic, referencePic, ... )
            %shift the DAT image:
            xShift = tform.T(3,1);
            yShift = tform.T(3,2);     
            [datTranslated,datTranslatedRef] = imtranslate(datIm,[xShift yShift],'OutputView','full');
            %TH pTranslated is the translated image; pTranslatedRef is an imref2D
            %object to tell MatLab where to place the image in a figure.
            RA = imref2d(size(aIm)); %TH creates the inref2D object for aMaxProj image. each of these objects are used by imshowpair, amung other fns, to place images correctly in figures
            %Plot the alignment, select the ANT and POST position in ref frame of zoomed out, anterior image:
            hAx = handles.axes1;
%             axis('off');
            imshowpair(datTranslated, datTranslatedRef,aIm,RA,'Scaling','independent','Parent',hAx);shg
            %add the aligned data image to align structure.
            align = get(handles.align,'UserData'); 
            align.datTranslated = datTranslated;
            align.datTranslatedRef = datTranslatedRef;
            %add the ANT to DAT FOR translational data to the align struc:
            align.xShift = -xShift; %negative values are bc the way I initially wrote it, we align the DATA im to the ANT im, so here we must switch the sign.
            align.yShift = -yShift;
         
            set(handles.align,'UserData',align) 
            toc
            
            set(handles.instrucs,'String','images aligned. once AP locations are saved you are ready to finish & save. Did alignment not work? try seeding it by selecting a corresponding point in Ant & Last T images.')
            set(handles.align,'String','Align')            
            
    end
elseif get(handles.seed,'Value') == 1 %use the pcik corresponding point feature (seeding)
    switch get(handles.alignMenu,'Value')
        case 1 %A-P
            tic;
            set(handles.instrucs,'String','aligning ANT & POST images. this may take a moment.')
            set(handles.align,'String','. . .')
            pause(0.1)
            %what channel:
            chanChoice = get(handles.chan1,'Value');
            if chanChoice == 1
                chan = 1;
            else 
                chan = 2;
            end
            %images:
            aIm = handles.antProjection(:,:,chan);
            pIm = handles.postProjection(:,:,chan);
            %rescale images for display:
            pIm = imadjust(mat2gray(pIm,[0 64000]));
            aIm = imadjust(mat2gray(aIm,[0 64000]));
            %first do initial hand alignment:
            %get the xs & ys
            antX = handles.antPoint(1);
            antY = handles.antPoint(2);
            postX = handles.postPoint(1);
            postY = handles.postPoint(2);
            %calc the hand shifts:
            xShiftHand = antX - postX;
            yShiftHand = antY - postY;
            %make crude alignment
            [pTranslatedHand,pTranslatedRefHand] = imtranslate(pIm,[xShiftHand yShiftHand],'OutputView','full');
            
            %TS, look at teh crude alignment:
%             hAx = handles.axes1;
%             RA = imref2d(size(aIm));
% %             imshow(imadjust(mat2gray(pTranslatedHand,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);
% %             axis('off');
%             imshowpair(pTranslatedHand, pTranslatedRefHand,aIm,RA,'Scaling','joint','Parent',hAx);shg
%             %store:
%             align.pTranslated = pTranslatedHand;
%             align.pTranslatedRef = pTranslatedRefHand;
%             align.aIm = aIm;
%             align.RA = RA;            
%             set(handles.align,'UserData',align) 
            
            %now do the acrtual alignment (this presently doesn't work. come back to it:
            %TH set up alignment params:
            [optimizer, metric] = imregconfig('multimodal');
            %here we reduce the InitialRadius property of the optimizer to ensure we
            %get an answer that converges
            r = optimizer.InitialRadius;
            optimizer.InitialRadius = r*10^(-3);
            optimizer.MaximumIterations = 500;
            optimizer.GrowthFactor = 1.5;
            optimizer.Epsilon = 1.5*10^(-10);
            %try messing with metric:
            metric.NumberOfSpatialSamples = 50000;
            metric.NumberOfHistogramBins = 500;
            % this gives you the transformation matrix
            RA = imref2d(size(aIm)); %TH creates the inref2D object for aMaxProj image. each of these objects are used by imshowpair, amung other fns, to place images correctly in figures
            tform = imregtform(pTranslatedHand, pTranslatedRefHand, aIm, RA, 'translation', optimizer, metric);       %syntax: imregtform(movingPic, referencePic, ... )
            %shift the POST image:
            xShift = tform.T(3,1);
            yShift = tform.T(3,2);     
            [pTranslated,pTranslatedRef] = imtranslate(pTranslatedHand,pTranslatedRefHand,[xShift yShift],'OutputView','full');
            %TH pTranslated is the translated image; pTranslatedRef is an imref2D
            %object to tell MatLab where to place the image in a figure.
            %Plot the alignment, select the ANT and POST position in ref frame of zoomed out, anterior image:
            hAx = handles.axes1;
            axis('off');
            imshowpair(pTranslated, pTranslatedRef,aIm,RA,'Scaling','independent','Parent',hAx);shg
            %set(handles.align,'UserData',[pTranslated, pTranslatedRef, aIm, RA]) 
            %store:
            align.pTranslated = pTranslated;
            align.pTranslatedRef = pTranslatedRef;
            align.aIm = aIm;
            align.RA = RA;            
            set(handles.align,'UserData',align) 
            
            toc
            set(handles.instrucs,'String','images aligned. now pick the ANT & POST positions. Note that this ')
            set(handles.align,'String','Align')

        case 2
            tic;
            set(handles.instrucs,'String','aligning last T point image to A-P image')
            set(handles.align,'String','. . .')
            pause(0.1)
            %what channel:
            chanChoice = get(handles.chan1,'Value');
            if chanChoice == 1
                chan = 1;
            else 
                chan = 2;
            end
            aIm = handles.antProjection(:,:,chan);
            datIm = handles.lastT; %this image onlyu has one channel
            %rescale images for display:
            datIm = imadjust(mat2gray(datIm,[0 64000]));
            aIm = imadjust(mat2gray(aIm,[0 64000]));
            %first do initial hand alignment:
            %get the xs & ys
            antX = handles.antPoint(1);
            antY = handles.antPoint(2);
            lastX = handles.lastTpoint(1);
            lastY = handles.lastTpoint(2);
            %calc the hand shifts:
            xShiftHand = antX - lastX;
            yShiftHand = antY - lastY;
            %make crude alignment
            [datTranslatedHand,datTranslatedRefHand] = imtranslate(datIm,[xShiftHand yShiftHand],'OutputView','full');  
            %TS:
            RA = imref2d(size(aIm));
            hAx = handles.axes1;
%             axis('off');
            imshowpair(datTranslatedHand, datTranslatedRefHand,aIm,RA,'Scaling','independent','Parent',hAx);shg
            %add the aligned data image to align structure.
            align = get(handles.align,'UserData'); 
            align.datTranslated = datTranslatedHand;
            align.datTranslatedRef = datTranslatedRefHand;
            %add the ANT to DAT FOR translational data to the align struc:
            align.xShift = -xShiftHand; %negative values are bc the way I initially wrote it, we align the DATA im to the ANT im, so here we must switch the sign.
            align.yShift = -yShiftHand;
         
            set(handles.align,'UserData',align)            
            
            %ATM this don't work too great. come back to it
%             %TH set up alignment params:
%             [optimizer, metric] = imregconfig('multimodal');
%             %here we reduce the InitialRadius property of the optimizer to ensure we
%             %get an answer that converges
%             r = optimizer.InitialRadius;
%             optimizer.InitialRadius = r*10^(-1);
%             optimizer.MaximumIterations = 500;
%             RA = imref2d(size(aIm)); 
%             % this gives you the transformation matrix
%             tform = imregtform(datTranslatedHand, datTranslatedRefHand, aIm, RS, 'translation', optimizer, metric);       %syntax: imregtform(movingPic, referencePic, ... )
%             %shift the POST image:
%             xShift = tform.T(3,1);
%             yShift = tform.T(3,2);     
%             [datTranslated,datTranslatedRef] = imtranslate(datTranslatedHand,datTranslatedRefHand,[xShift yShift],'OutputView','full');
%             %TH pTranslated is the translated image; pTranslatedRef is an imref2D
%             %object to tell MatLab where to place the image in a figure.
%             %Plot the alignment, select the ANT and POST position in ref frame of zoomed out, anterior image:
%             hAx = handles.axes1;
%             axis('off');
%             imshowpair(datTranslated, datTranslatedRef,aIm,RA,'Scaling','joint','Parent',hAx);shg
%             %add the aligned data image to align structure.
%             align = get(handles.align,'UserData'); 
%             align.datTranslated = datTranslated;
%             align.datTranslatedRef = datTranslatedRef;
%             %add the ANT to DAT FOR translational data to the align struc:
%             align.xShift = -xShiftHand; %negative values are bc the way I initially wrote it, we align the DATA im to the ANT im, so here we must switch the sign.
%             align.yShift = -yShiftHand;
%          
%             set(handles.align,'UserData',align) 
            toc
            
            set(handles.instrucs,'String','images aligned. once AP locations are saved you are ready to finish & save.')
            set(handles.align,'String','Align')            
            
    end
    
elseif get(handles.manual,'Value') == 1
    switch get(handles.alignMenu,'Value')
        case 1 %A-P
            tic;
            set(handles.instrucs,'String','aligning ANT & POST images')
            set(handles.align,'String','. . .')
            pause(0.1)
            %what channel:
            chanChoice = get(handles.chan1,'Value');
            if chanChoice == 1
                chan = 1;
            else 
                chan = 2;
            end
            %images:
            aIm = handles.antProjection(:,:,chan);
            pIm = handles.postProjection(:,:,chan);
            %rescale images for display:
            pIm = imadjust(mat2gray(pIm,[0 64000]));
            aIm = imadjust(mat2gray(aIm,[0 64000]));
            %first do initial hand alignment:
            %get the xs & ys
            antX = handles.antPoint(1);
            antY = handles.antPoint(2);
            postX = handles.postPoint(1);
            postY = handles.postPoint(2);
            %calc the hand shifts:
            xShiftHand = antX - postX;
            yShiftHand = antY - postY;
            %make crude alignment
            [pTranslatedHand,pTranslatedRefHand] = imtranslate(pIm,[xShiftHand yShiftHand],'OutputView','full');
            %look at teh alignment:
            hAx = handles.axes1;
            RA = imref2d(size(aIm));
            imshowpair(pTranslatedHand, pTranslatedRefHand,aIm,RA,'Scaling','independent','Parent',hAx);shg
            %store:
            align.pTranslated = pTranslatedHand;
            align.pTranslatedRef = pTranslatedRefHand;
            align.aIm = aIm;
            align.RA = RA;            
            set(handles.align,'UserData',align) 

            toc
            set(handles.instrucs,'String','images aligned. now pick the ANT & POST positions.')
            set(handles.align,'String','Align')

        case 2
            tic;
            set(handles.instrucs,'String','aligning last T point image to A-P image')
            set(handles.align,'String','. . .')
            pause(0.1)
            %what channel:
            chanChoice = get(handles.chan1,'Value');
            if chanChoice == 1
                chan = 1;
            else 
                chan = 2;
            end
            aIm = handles.antProjection(:,:,chan);
            datIm = handles.lastT; %this image onlyu has one channel
            %rescale images for display:
            datIm = imadjust(mat2gray(datIm,[0 64000]));
            aIm = imadjust(mat2gray(aIm,[0 64000]));
            %first do initial hand alignment:
            %get the xs & ys
            antX = handles.antPoint(1);
            antY = handles.antPoint(2);
            lastX = handles.lastTpoint(1);
            lastY = handles.lastTpoint(2);
            %calc the hand shifts:
            xShiftHand = antX - lastX;
            yShiftHand = antY - lastY;
            %make crude alignment
            [datTranslatedHand,datTranslatedRefHand] = imtranslate(datIm,[xShiftHand yShiftHand],'OutputView','full');  
            %TS:
            RA = imref2d(size(aIm));
            hAx = handles.axes1;
%             axis('off');
            imshowpair(datTranslatedHand, datTranslatedRefHand,aIm,RA,'Scaling','independent','Parent',hAx);shg
            %add the aligned data image to align structure.
            align = get(handles.align,'UserData'); 
            align.datTranslated = datTranslatedHand;
            align.datTranslatedRef = datTranslatedRefHand;
            %add the ANT to DAT FOR translational data to the align struc:
            align.xShift = -xShiftHand; %negative values are bc the way I initially wrote it, we align the DATA im to the ANT im, so here we must switch the sign.
            align.yShift = -yShiftHand;
         
            set(handles.align,'UserData',align)            
            
            toc
            
            set(handles.instrucs,'String','images aligned. once AP locations are saved you are ready to finish & save.')
            set(handles.align,'String','Align')            
            
    end
        
end
        


% --- Executes on button press in pickApoint.
function pickApoint_Callback(hObject, eventdata, handles)
% hObject    handle to pickApoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(handles.displayMenu,'Value')
    case 1
        [xpt ypt]=ginput(1);
        handles.antPoint = [xpt ypt];
        
        set(handles.instrucs,'String','ANT point saved')
    case 2
        [xpt ypt]=ginput(1);
        handles.postPoint = [xpt ypt];
        
        set(handles.instrucs,'String','POST point saved')
    case 3
        [xpt ypt]=ginput(1);
        handles.lastTpoint = [xpt ypt];
        
        set(handles.instrucs,'String','Last T point saved')
    case 4
        set(handles.instrucs,'String','invalid menu selection, pick ANT, POST or LAST T.')
        return
    case 5
        set(handles.instrucs,'String','invalid menu selection, pick ANT, POST or LAST T.')
        return
end
guidata(hObject, handles);

% --- Executes on button press in seed.
function seed_Callback(hObject, eventdata, handles)
% hObject    handle to seed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of seed


% --- Executes on button press in pickAP.
function pickAP_Callback(hObject, eventdata, handles)
% hObject    handle to pickAP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.instrucs,'String','select first the ANT, then the POST, right click when done')
pause(0.1)

hold off
%first show the aligned images to clear any previous AP selections:
align  = get(handles.align,'UserData');
if isempty(align)
    set(handles.instrucs,'String','no aligment info. Make one first.')
    return
end
pTranslated = align.pTranslated;
pTranslatedRef = align.pTranslatedRef;
aIm = align.aIm;
RA = align.RA;
%show it:
hAx = handles.axes1;
imshowpair(pTranslated, pTranslatedRef, aIm, RA,'Scaling','joint','Parent',hAx);shg
%select ANT & POST positions, with some bells and whistles:
axes(handles.axes1)
[xs, ys, but]=ginputc(3,'Color','b','ShowPoints',true, 'ConnectPoints', true);
set(handles.figure1,'waitstatus','waiting') %for passing AP coords to flimscroll
axes(handles.axes1)
hold on;
% plot(xs(1:2), ys(1:2), '-+r','LineWidth',1,'MarkerSize',12,'Parent',hAx);shg
%labeling the selections:
txt1 = ['ANT'];
txt2 = ['POST'];
text(xs(1),ys(1),txt1);
text(xs(2),ys(2),txt2);

plot(xs(1:2), ys(1:2), '-+r','LineWidth',1,'MarkerSize',12);shg %,'Parent',hAx
%AP positons in the ref frame of the stiched together ANT & POST images:
AP.coordA = [xs(1) ys(1)];
AP.coordP = [xs(2) ys(2)];

%store:
set(handles.pickAP,'UserData',AP)
set(handles.instrucs,'String','AP position saved')


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%this is to get the GUI to output (and pass) the coord struc to flimscroll
if isequal(get(hObject,'waitstatus'), 'waiting')
    uiresume(hObject);
else
    delete(hObject);
end
%default:
% delete(hObject)

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get AP in ANT ref frame:
antAPcoords = get(handles.pickAP,'UserData');
if isempty(antAPcoords)
    set(handles.instrucs,'String','no AP coordinates, pick those then come back to this.')
    return
end
antA = antAPcoords.coordA; %[xpos ypos]
antP = antAPcoords.coordP;
%get the translation info from the LAST T to ANT alignment:
align = get(handles.align,'UserData');
if isempty(align)
    set(handles.instrucs,'String','no last T point to AP alignment data, get that then come back to this.')
    return
end
xShift = align.xShift; %a scaler
yShift = align.yShift;
%Get the AP position in the LAST T FOR:
%have to swict the y axis sign bc of the weird way matlab plots images,
%with the origin in the upper LH corner.
coordA = [antA(1)+xShift antA(2)+yShift];
coordP = [antP(1)+xShift antP(2)+yShift];
%put these into a struc for output
coordStruc.coordA = coordA;
coordStruc.coordP = coordP;
%save it to GUI data output:
set(handles.save,'UserData',coordStruc);
%save coords to a file:
currDir = cd;
cd(handles.dir{1})
% cd ..
fp = cd;
fn = extractBefore(handles.dir{2},'_'); %try this, should work well enough to get some info on the image to the filename    
fileOutPath = [fp fn '-APpositionStructure.dat'];
eval(['save ' fileOutPath ' coordStruc']);
cd(currDir)

% %TS, show the AP plotted on the data image to confirm
% hold off
% Imax = handles.lastT;
% hAx = handles.axes2;
% axes(hAx);
% plot(coordA(1),coordA(2),'bx');
% hold on
% plot(coordP(1),coordP(2),'rx');shg
% imshow(imadjust(mat2gray(Imax,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);shg
% plot(coordA(1),coordA(2),'yx');
% plot(coordP(1),coordP(2),'rx');shg
% axis('on')

set(handles.instrucs,'String','all good to go. AP location saved in ...-APlocationStructure.dat file. close this window to pass AP location & return to flimscroll.')


% --- Executes on button press in display2.
function display2_Callback(hObject, eventdata, handles)
% hObject    handle to display2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of display2


% --- Executes on button press in manual.
function manual_Callback(hObject, eventdata, handles)
% hObject    handle to manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of manual


% --- Executes on selection change in alignMenu.
function alignMenu_Callback(hObject, eventdata, handles)
% hObject    handle to alignMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns alignMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from alignMenu


% --- Executes during object creation, after setting all properties.
function alignMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alignMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in displayIn2.
function displayIn2_Callback(hObject, eventdata, handles)
% hObject    handle to displayIn2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hold off
%display in axes 2:
hAx = handles.axes2;

switch get(handles.displayMenu,'Value')
    case 1 %anterior
        Imax = handles.antProjection;
        if isempty(Imax)
            set(handles.instrucs,'String','no max proj image. Make or load one')
            return
        end
        %show selected channel:
        chanChoice = get(handles.chan1,'Value');
        if chanChoice == 1
            chan = 1;
        else
            chan = 2;
        end
        im = Imax(:,:,chan);
        clowval = min(min(im));
        chival = max(max(im));
        axes(hAx);
        imagesc(im,[clowval chival]);axis(hAx,'equal','off');colormap(gray(256));
%         imshow(imadjust(mat2gray(im,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);
    case 2 %posterior
        Imax = handles.postProjection;
        if isempty(Imax)
            set(handles.instrucs,'String','no max proj image. Make or load one')
            return
        end

        chanChoice = get(handles.chan1,'Value');
        if chanChoice == 1
            chan = 1;
        else
            chan = 2;
        end
        im = Imax(:,:,chan);
        clowval = min(min(im));
        chival = max(max(im));
        axes(hAx);
        imagesc(im,[clowval chival]);axis(hAx,'equal','off');colormap(gray(256));
%         imshow(imadjust(mat2gray(im,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);
    case 3 %last T
        Imax = handles.lastT;
        if isempty(Imax)
            set(handles.instrucs,'String','no max proj image. Make or load one')
            return
        end
        clowval = min(min(Imax));
        chival = max(max(Imax));
        axes(hAx);
        imagesc(Imax,[clowval chival]);axis(hAx,'equal','off');colormap(gray(256));
%         imshow(imadjust(mat2gray(Imax,[0 64000])),'DisplayRange',[],'InitialMagnification',100,'Parent',hAx);
        
    case 4 %A-P alignment
        %get the aligned images:
        align  = get(handles.align,'UserData');
        if isempty(align)
            set(handles.instrucs,'String','no aligment info. Make one first.')
            return
        end
        pTranslated = align.pTranslated;
        pTranslatedRef = align.pTranslatedRef;
        aIm = align.aIm;
        RA = align.RA;
        %show it:
        %hAx = handles.axes1;
        imshowpair(pTranslated, pTranslatedRef, aIm, RA,'Scaling','joint','Parent',hAx);shg
%         axis('off');
        %add AP axis if it exists:
        AP = get(handles.pickAP,'UserData');
        if ~isempty(AP)
            hold on;
            xs = [AP.coordA(1) AP.coordP(1)];
            ys = [AP.coordA(2) AP.coordP(2)];
            %labeling the selections:
            txt1 = ['ANT'];
            txt2 = ['POST'];
            text(xs(1),ys(1),txt1);
            text(xs(2),ys(2),txt2);
            plot(xs(1:2), ys(1:2), '-+r','LineWidth',1,'MarkerSize',12);shg %,'Parent',hAx
        end
        
        set(handles.instrucs,'String','AP alignment displayed')
        
    case 5 %AP-data alignment
        %get the aligned images:
        align  = get(handles.align,'UserData');        
        if isempty(align)
            set(handles.instrucs,'String','no aligment info. Make one first.')
            return
        end
        datTranslated = align.datTranslated;
        datTranslatedRef = align.datTranslatedRef;
        aIm = align.aIm;
        RA = align.RA;
        %show it:
        %hAx = handles.axes1;
        imshowpair(datTranslated, datTranslatedRef, aIm, RA,'Scaling','joint','Parent',hAx);shg
%         axis('off');
        set(handles.instrucs,'String','AP-data alignment displayed')        
        
end
