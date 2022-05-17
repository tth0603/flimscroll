function mih = firstEventBootstrapGenerator2(varargin)
%   mih = firstEventBootstrapGenerator2(analData1,analData2,...,bins,binsOfInterest,numOfBS)
%
% outputs the bootstrapped data for use with firstEventBootstrapPlotter.m &
% firstEventBootstrapFitter.m. does not shift data along x axis. 
% 
% here we grab a nuke, not a first passage time, to get the estimate for
% the active fraction as well as the rate/steps
%
% INPUTS
%   cia == cumulative interval array from flimscroll2
%   numOfBS == number of bootstrap data replicates
%
% OUTPUTS
%   mib == an numOfBS x numOfEvents matrix containing the bootstrapped
%   data. you can use this mat to fit for fit param error and plot for
%   curve error bars
%
% Timothy Harden 2020

% initialize the data:
numberOfReps = length(varargin) - 3;

%initialize some variables:
bins = varargin{numberOfReps + 1};
boi = varargin{numberOfReps + 2};
numOfBS = varargin{numberOfReps + 3};

% get the matrix that we'll be pulling boot straps from. This is a cia from
% both reps with nukes in the bins of interest.

firstPassageV = []; % a list of all first passage times within our bins of interest for both replicates. the zeros represent nukes within the bins that never saw an event
for i = 1:numberOfReps
    % binning shit:
    out = positionDistribution(varargin{i},bins);
    % get first passage times in cia form:
    cia = varargin{i}.cia;
    firstMat = cia(cia(:,2) == -3,:);
    % Nuke position mat [1.nukeNumber 2.nukePostion 3.binAssignment]
    eventPos = out.eventPos;
    % Nukes without an event mat [1.nukeNumber 2.nukePostion 3.binAssignment]
    noEventNukePos = out.noEventNukePos;
    % get the nukes in the bins of interest
    nukes = [];
    noEventNukes = [];
    for j = 1:length(boi)
        nukes = [nukes;eventPos(eventPos(:,3) == boi(j),1)];
        noEventNukes = [noEventNukes;noEventNukePos(noEventNukePos(:,3) == boi(j),1)];
    end
    % get the first passage times for these nukes:
    for j = 1:length(nukes)
        firstPassageV = [firstPassageV;firstMat(firstMat(:,1) == nukes(j),7)];
    end
    % tack zeros for the nukes that saw no event:
    firstPassageV = [firstPassageV;noEventNukes*0];
end

% BS part. you have a list of first passage events (and zeros) that you're BSing. The length of this list
% is what you're normalizing to. 

%get the BS mat
intervals = firstPassageV;
n = size(firstPassageV,1);
%make the BS mat - cols are BS, rose are time points
tbl=[[1:n]' 1/n*ones(n,1)];  %this is from LJF's method of picking with replacement
bsMat = zeros(n,numOfBS);
for i = 1:numOfBS
    indx = probability_steps(tbl,n);
    bs = intervals(indx); 
    bsMat(:,i) = bs;
end

mih = bsMat;