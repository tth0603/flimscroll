function mih = firstEventBootstrap2(varargin)
%   [mib,plotHandle] = firstEventBootstrap(cia,normalization,numOfBS,plotColor,fignum)
% My latest and greatest bootstrap script. will spit out and plot the 5th &
% 95th percentile bootstrap vector so as to get error bars.
%
% INPUTS
%   cia == cumulative interval array from flimscroll2
%   normalization == scaler. what to normalize to. can be empty, and the curve is not normalized
%   numOfBS == number of bootstrap data replicates
%
% OUTPUTS
%   mib == an numOfBS x numOfEvents matrix containing the bootstrapped data
%
% Timothy Harden 2020

%initiate:
cia = varargin{1};
normalization = varargin{2};
numOfBS = varargin{3};

%grab the first event (its length)
first = cia(cia(:,2) == -3,7); %a vector
%make the first event occur at t = 0;
first = (first - min(first) + 1); %harden 202009: you can't have your first event start at zero. or maybe you can. whatever. 
%determine normalization
if ~isempty(normalization)
    norm = normalization;
else 
%     mx = length(first); %there is a debate in my head as to which is more
%     useful
    norm = 1;
end

%get the BS mat
intervals = first;
n = size(first,1);
mx = ceil(max(intervals)); %the longest event in data
%make the BS mat - cols are BS, rose are time points
tbl=[[1:n]' 1/n*ones(n,1)];  %this is from LJF's method of picking with replacement
bsMat = zeros(n,numOfBS);
parfor i = 1:numOfBS
    indx = probability_steps(tbl,n);
    bs = intervals(indx); 
    bsMat(:,i) = bs;
end

mih = bsMat;