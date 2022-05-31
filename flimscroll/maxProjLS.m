function maxProjLS(acqnDate,dataDir,dirVect,varargin)
% maxProjLS(acqnDate,dataDir,dirVect,alternating,fileName,gfpChan)
%
% this fn will max proj your LS data images, then you can use flimscroll2 to
% process & analyze them
%
% This script should come with the file filelocations.dat.
% This script will go directory by direcory writing the MS2 max proj image
% files. If these files already exist it will skip writing them. It took 
% 8.5 seconds to run this script, pointing at two directories, to get to 
% ms2scroll without writing any images. It will then call a GUI (ms2scroll) 
% to define a ROI that will be analyzed on down the pipeline
%
% for 216 time points this took ~1hr to run writing images with the
% original files on TH D:\ drive
%
% Start by editing filelocations.dat:
% 1. load it (run: load C:\Users\tth12\matlab\fig_files\flimscroll\fliFileLocations.dat -mat) 
% (change the file path)
% 2. This contains a structure (fileLocations) with members (each a string):
%       a. inPath: directory where your data files are stored
%       b. outPath: where you'd like to write the output
% 3. edit each of these. you should only have to do this once if you keep
% your directories consistent.
%   expl: 
%           fileLocations.inPath = 'D:\matlab\image-data\'
%           fileLocations.outPath = 'C:\Users\tth12\matlab\data\LSanal\'
%       Then save the new locations to filelocations.dat:
%           save C:\Users\tth12\matlab\fig_files\LSdataAnal\fileLocation.dat fileLocations
%
%   INPUTS (these are the same as for exportDataForFishTH):
%       acqnDate == directory data is containted in. expl: date = '2016-10-20'
%       dataDir == subdirectory within 'acqnDate' containing specific acqn. expl: dataDir = 'zelda1'
%       dirVect == a vector specifying the numbers of the subdirectories
%       you want to process in the dataDir directory. Exp: [1:5] or [2:3]
%       or [1,3:5]. you can pass this argument with '[]'. 
%       
%       optional inputs:
%           alternating == if the acqn is fucked and the channels alternate
%           between times, specify here with the string 'alt'. If also
%           specifying gfpChan, spcify which is the gfpChan of the first
%           time point
%
%           fileName == as for gfpChan, but a cell array of strings of
%           fileNames. expl: fileName{1} = 'MMStack_Pos0.ome.tif'; etc...
%           or:
%           fileName = ['aqn1_MMStack_Pos0.ome.tif' 'aqn1_MMStack_Pos0.ome.tif' 'aqn3_MMStack_Pos0.ome.tif'];
%           fileName = cellstr(fileName);
%
%           gfpChan == a vector of the gfp channels for each of the acqns
%           indicated in dirVect. Must be same length as dirVect.expl: [1 1
%           1 1 1] or [1 2] or [1 2 2 2]. Presently, the script figures out
%           the gfp channel for itself, so this input is unnessary in all
%           cases
%
%   USAGE:  maxProjLS(acqnDate,dataDir,[1:5],fileName,'alt',[1 1 1 1 1])
%           or
%           maxProjLS(acqnDate,dataDir,[1:5])
%           or 
%           maxProjLS(acqnDate,dataDir,[])
%           or
%           maxProjLS(acqnDate,dataDir,[],'alt')
%
%   OUTPUT:
%       none at present.
%
% Timothy Harden 2018
%
%for future: 
% 1. make a test case, so you dont have to waste too much time
%TSing. and along this line, it may not be best to stick all tifs in one
%file... There is code to write each max proj to its own file.
% 2. fix the clear heap space protocol. Matlab 2018a allows heap space up
% to 15GB, which works well for my machine. Other machines may struggle and
% older versions of matlab only allow 8GB. See https://www.mathworks.com/matlabcentral/fileexchange/36757-java-heap-cleaner

%% some ome retrieval stuff:
%trying to get channel colors: 
% %green?:
% chan1col = omeMeta.getChannelColor(0,0).getValue(); %indicies: (imageIndex, chanIndex)
% %red?:
% chan2col = omeMeta.getChannelColor(0,1).getValue();
% chanID = omeMeta.getChannelID(0,1); %same indicies as above

% the channel for a given plane:
% planeChan = omeMeta.getPlaneTheC(0,798); %(imageIndex,planeIndex); base zero
% Time
%  planeChan = omeMeta.getPlaneTheT(0,798); %(imageIndex,planeIndex); base zero
% number of planes:
% planeCount = omeMeta.	getPlaneCount(0); %(imageIndex)

%% preamble from export data for fish TH
tic;
initialTimePoint = 1;   %may have to make this an input 
if isempty(dirVect)
    dirVect = [1];
end
len = length(dirVect);
%initiate variables:
alt = false;
for j = 1:len
    fileName{j} = 'MMStack_Pos0.ome.tif';
end
gfpChan = ones(len,1);
findGfp = true;
%sort optional inputs:
if nargin > 3
    ln = length(varargin);
    for i = 1:ln
        if strcmp(varargin{i},'alt')
            alt = true;
        end
        if iscell(varargin{i})
            clear fileName
            fileName = varargin{i};
        end 
        if ~ischar(varargin{i}) && ~iscell(varargin{i})
            gfpChan = varargin{i};
            findGfp = false;
        end
    end
end

%make the his chan the other one:
for i = 1:len
    if gfpChan(i) == 1
        hisChan(i) = 2;
    elseif gfpChan(i) == 2
        hisChan(i) = 1;
    else
        warning('The channels are funny. are there more than 2?')
    end
end

load fliFileLocations.dat -mat % a pre defined file. see help above

%if the ouput directories don't exist, make them:
maxPath = [fileLocations.outPath acqnDate filesep dataDir];
%PreProcPath = [fileLocations.outPath 'PreProcessedData' filesep 'test_5' filesep acqnDate filesep dataDir];
if ~exist(maxPath,'dir')
    mkdir(maxPath)
end

%do we want to make ini time point an input?
if isempty(initialTimePoint)
    initialTimePoint = 1;
end

%for clearing heap space (needed for jheapcl.m) ... doesn't seem to work
% jp = fileLocations.javaPath;
% javaaddpath(which(jp))

for i = 1:len
    %define dir
    fileInPath = [fileLocations.inPath acqnDate filesep dataDir filesep fileName{i}];
    fileInPath = char(fileInPath);
    %define nT:
    %have to bfopen first fo nT:
    reader = bfGetReader(fileInPath); %low level retireival fn to get file reader without loading the data
    omeMeta = reader.getMetadataStore(); %retrieves ome metadata
    nT = omeMeta.getPixelsSizeT(0).getValue(); %number of time points
    nZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
    nX = omeMeta.getPixelsSizeX(0).getValue();
    nY = omeMeta.getPixelsSizeY(0).getValue();
    nC = omeMeta.getPixelsSizeC(0).getValue();
    %get the gfpChan (IN BETA)
    if findGfp
        chanVal = omeMeta.getChannelColor(0,0).getValue();
        if chanVal < -5000000 %empirically determined. This may vary and thereby need to be changed. 
            gfpChan(i) = 1; %chan 1 is grn
            hisChan(i) = 2;
        else
            gfpChan(i) = 2; %chan 2 is grn
            hisChan(i) = 1;
        end
    end            
    %preallocate: (this loop took 32s for 6 images
    Imax = ones(nY,nX,nT);
    %Imax = [];  %TSing
    fileOutPath = [maxPath filesep acqnDate '-' dataDir '-MS2proj.tiff'];  %this slips the name of the input folder in front of each file name. depending on what you want to do with the files next (imaris, HG pipeline, etc) this may need to be revised
    %fileOutPath = [maxPath filesep acqnDate '-' dataDir '-MS2projLastTest.tiff'];    %TSing
    for iT = 1:nT %one iteration takes about 4s
    %for iT = 1:23    %TSing
        fprintf('%d MS2 projection \n',iT + initialTimePoint - 1)
        %for alternating channels, switch the channel for even time points0
        if alt && mod(iT,2) == 0 %mod thingy is zero for even ints
            chan = hisChan(i);
        else 
            chan = gfpChan(i);
        end
        I = ones(nY,nX,nZ);    
        %I = [];   %TSing
        for iZ = 1:nZ
            iPlane = reader.getIndex(iZ - 1, chan - 1, iT - 1) + 1; 
            I(:,:,iZ) = bfGetPlane(reader, iPlane);
        end
        Imax(:,:,iT) = max(I,[],3);
        fprintf('.')
        fprintf('\n')
%         %% TSing
%         if iT == 21
%             test = 1;
%         end
%         %% TSing
%         iT = 21;
%         for iZ = 1:nZ
%             iPlane = reader.getIndex(iZ - 1, chan - 1 , iT - 1) + 1; 
%             I(:,:,iZ) = bfGetPlane(reader, iPlane);
%         end
%         Imax(:,:,iT) = max(I,[],3);
%         figure(32);imshow(imadjust(mat2gray(Imax(:,:,iT),[0 64000])),'DisplayRange',[],'InitialMagnification',100);
%         %% 
%         iPlane = reader.getIndex(0, 0, 20) 
%         tst = bfGetPlane(reader, 3201);
%         figure(40);imshow(imadjust(mat2gray(tst,[0 64000])),'DisplayRange',[],'InitialMagnification',100);
%         %%
% %         reader2 = bfGetReader(fileInPath);
%         iPlane = reader2.getIndex(40, 0, 19) + 1
%         tst = bfGetPlane(reader2, 3201);
%         figure(40);imshow(imadjust(mat2gray(tst,[0 64000])),'DisplayRange',[],'InitialMagnification',100);
%         %%
        %%%%%clearing java heap space:
        %in previous runs, I discovered for typical LS imaging conditions,
        %we run out of heap space (set at 8.168 GB)
%         javaClear = iT/50; %clear java heap space every 50 frames
%         if isinteger(javaClear)
%             jheapcl
%         end
        %%%%%%%%%%%doesn't seem to work
    end
%for multiple acquisitions, this breaks the script. Not sure how much I
%care abotu this now:
    bfsave(Imax,fileOutPath, 'BigTiff', true); %need big tiff! 

clear I 
    %max proj the his images:
    hisImax = ones(nY,nX,nT);
    %hisMax = [];
    hisOutPath = [maxPath filesep acqnDate '-' dataDir '-HISproj.tiff'];  
    %hisOutPath = [maxPath filesep acqnDate '-' dataDir '-HISprojLastTest.tiff'];  %use this for TSing
    for iT = 1:nT %one iteration takes about 4s
    %for iT = nT   %TSing
        fprintf('%d HIS projection \n',iT + initialTimePoint - 1)
        %for alternating channels, switch the channel for even time points0
        if alt && mod(iT,2) == 0 %mod thingy is zero for even ints
            chan = gfpChan(i);
        else 
            chan = hisChan(i);
        end
        I = ones(nY,nX,nZ); 
        %I = [];  %TSing
        for iZ = 1:nZ
            iPlane = reader.getIndex(iZ - 1, chan - 1, iT - 1) + 1; 
            I(:,:,iZ) = bfGetPlane(reader, iPlane);
        end
        hisImax(:,:,iT) = max(I,[],3);
        fprintf('.')
        fprintf('\n')
%         javaClear = iT/50; %clear java heap space every 50 frames
%         if isinteger(javaClear)
%             jheapcl
%         end
    end
    bfsave(hisImax,hisOutPath, 'BigTiff', true);

    initialTimePoint = initialTimePoint + nT;  %make sure this iterative indexing is correct ...it is
end


toc
