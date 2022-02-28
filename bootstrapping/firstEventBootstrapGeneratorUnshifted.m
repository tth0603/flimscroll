function mih = firstEventBootstrapGeneratorUnshifted(varargin)
%   mih = firstEventBootstrapGenerator(cia,numOfBS)
%
% outputs the bootstrapped data for use with firstEventBootstrapPlotter.m &
% firstEventBootstrapFitter.m
%
% INPUTS
%   cia == cumulative interval array from flimscroll2
%   numOfBS == number of bootstrap data replicates
%
% OUTPUTS
%   mib == an numOfBS x numOfEvents matrix containing the bootstrapped data
%
% Timothy Harden 2020

%initiate:
cia = varargin{1};
numOfBS = varargin{2};

%grab the first event (its length)
first = cia(cia(:,2) == -3,7); %a vector

%get the BS mat
intervals = first;
n = size(first,1);
mx = ceil(max(intervals)); %the longest event in data
%make the BS mat - cols are BS, rose are time points
tbl=[[1:n]' 1/n*ones(n,1)];  %this is from LJF's method of picking with replacement
bsMat = zeros(n,numOfBS);
for i = 1:numOfBS
    indx = probability_steps(tbl,n);
    bs = intervals(indx); 
    bsMat(:,i) = bs;
end

mih = bsMat;