function timebaseStruc = makeTimeBaseLS(varargin)
%
% timebaseStruc = makeTimeBaseLS(numberOfAcqusitions,save)
%
% This script can take a while. A couple mins to 10s of mins if your data
% is big or far away or your computer sucks.
%
% This will allow you to make a time base from one or a series of diSPIM
% acquisitions. This can be run alongside MaxProjLS prior to flimscroll or
% called from within flimscroll via Make time base. The latter will only 
% handle single acquisitions and must be saved from flimscroll.
% The user will be prompted to select the first image file for each 
% acquisition to be included in the time base. This script calls a 
% fileLocations structure with the fliFileLocations.dat file contained in the 
% Matlab file path.
%
% INPUTS:
%       numberOfAcqusitions == Number of image acqns to include. An
%       optional scaler
%       save == an optional string specifying if you'd like to save the
%       time base to a .dat file. (e.g. 'save' , 'S' , 'any string really')
%       These inputs must be put in order if both are included.
%
% OUTPUTS:
%       mib = a structure with two members: the time base itself and a
%       vector of the delta T's between acqns. 
%
%       The script will also save the output structure to a .dat file as
%       directed but the outPath member of the fileLocations struc.
%
% USAGE:
%       timebaseStruc = makeTimeBaseLS(2,'save');
%       timebaseStruc = makeTimeBaseLS;
%
% Harden 2018

%have user pick all the acquisitions;
switch length(varargin)
    case 0
        len = 1;
        sv = false;
    case 1
        if isscalar(varargin{1})
            len = varargin{1}; 
            sv = false;
        elseif ischar(varargin{1})
            len = 1;
            sv = true;
        else 
            disp('inputs are funny. try again.')
            return
        end
    case 2
        len = varargin{1};
        sv = true;
end

if sv
    try
        eval(['load fliFileLocations.dat -mat'])
    catch
        fprintf('no fliFileLocations.dat file. add that to the file path and come back \n')
        return
    end
end

currDir = cd;
for i = 1:len
    [fn fp] = uigetfile('*.*','pick the first file in the acqn series');
    readerCell{i} = bfGetReader([fp fn]);
    cd(fp)
    cd ..
end
cd(currDir)

%new for loop to not keep the user waiting forever:
%now build the time base, based on ometoimarisTTB.m
iniT = [];
deltaV = [];
ttb = [];
for i = 1:len
    %get the bioformats stuff
    reader = readerCell{i};
    omeMeta = reader.getMetadataStore();
    
    nZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
    nC = omeMeta.getPixelsSizeC(0).getValue(); % number of channels
    nT = omeMeta.getPixelsSizeT(0).getValue(); % number of time points
    
    %for time base:
    dateStr = omeMeta.getImageAcquisitionDate(0).getValue();  %TH get the initial time af acqn
    dateT = datetime(dateStr,'InputFormat','yyyy-MM-dd''T''HH:mm:ss');  %TH make it a datetime matlab thingy
    dateV = datevec(dateT);     %TH convert it to a time vecctor
    iniT = [iniT; dateV];        %TH make a vector the absolute initial times
    deltaT = etime(dateV, iniT(1,:)); %TH find the diff between the first acqn and current acqn
    deltaV = [deltaV; deltaT];      %TH make a vector of the times from the initial acqn (This should be one of the outputs of the script)
    
    dT = [];    
    for j = 0:nC*nZ:nC*nT*nZ-1 %TH this plucks the time of the first image of each stack for channel one at each time point. most omeMeta stuff starts indexing at 0
        planeT = double(value(omeMeta.getPlaneDeltaT(0,j)));  %TH for loops makes your time vector in ms
        dT = [dT; planeT];
    end
    
    dTsec = dT./1000; %TH convert time to seconds (what etime outputs)
    relT = dTsec + deltaT; %TH add the time between acqns
    ttb = [ttb; relT];
    
    % iSeries = 1;
    % reader.setSeries(iSeries - 1);

end
%output it
timebaseStruc.ttb = ttb;
timebaseStruc.deltaV = deltaV;

%save it:
if sv 
    outPath = fileLocations.outPath;
    newFn = extractBefore(fn,'_');
    eval(['save ' outPath newFn '-ttb.dat timebaseStruc']);
    fprintf('time base saved to %s%s-ttb.dat \n',outPath,newFn)
end


