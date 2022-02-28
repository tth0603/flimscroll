function [mih,plotHandle] = firstEventBootstrapPlotter(varargin)
%   [mih,plotHandle] = firstEventBootstrapPlotter(bsMat,normalization,plotColor,fignum)
%
% Spits out and plots the 5th & 95th percentile bootstrap vector so as to get error bars.
%
% INPUTS
%   BSmat == matstrix of bootstrapped data from firstEventBootstrapGenerator
%   norm == a scaler. what we normalize to. Ignores if empty ([])
%   colr == vetor length 3
%   figNum == scaler
%
% OUTPUTS
%   mib == an N x 2 matrix of the upper and lower sints, respectively
%   plotHandle == plot handle for editing the ploted curves post hoc
%   It will also plot a shaded area replresenting the error "bars" on a
%   natural log scale
%
% Timothy Harden 2020

%initiate:
bsMat = varargin{1};
normalization = varargin{2};
colr = varargin{3};
figNum = varargin{4};

mx = ceil(max(max(bsMat)));  %mx = ceil(max(intervals)); %the longest event in data
numOfBS = size(bsMat,2);
%determine normalization
if ~isempty(normalization)
    norm = normalization;
else 
%     mx = length(first); %there is a debate in my head as to which is more
%     useful
    norm = 1;
end


%sort the individual bootstraps (ie the cols of bsMat)
sintsMat = zeros(mx + 1,numOfBS); %plus one bc we start from zero time
for j = 1:numOfBS
    sints = [];
    for i = 0:mx %zero bc time
        logik = bsMat(:,j) < i;
        sints = [sints; sum(logik)];
    end
    sintsMat(:,j) = sints;
end
%normailization(?)
sintsMat = sintsMat./norm;

%so this is a bit tricky to wrap my mind around, but i'll try:
% for each time (row) we will sort all the different BS vector values in
% ascending order. then we can pluckl out the 5th & 95th percentile at each
% time. make sense?
B = sort(sintsMat,2);

%pluck out the 5th & 95th:
lowerI = floor(0.05*numOfBS); %5th index
if lowerI == 0
    lowerI = 1;
end
upperI = ceil(0.95*numOfBS);
lower = B(:,lowerI);
upper = B(:,upperI);

%time vector:
t = [0:mx]';

bnd = size(upper,1);
lower = lower';
upper = upper'; %bc they need to be row v's for some reason
ts = [0:bnd - 1, fliplr(0:bnd - 1)]; %not sure why we gotta flip these things or if we really need fliplr and not just flip but whatever
% shade = [log(lower2/n), fliplr(log(upper2/n))];
% shade = [lower2,fliplr(upper2)];
shade = [lower,fliplr(upper)];
figure(figNum);
plotHandle = fill(ts,shade,colr); 
set(plotHandle,'facealpha',0.25,'edgealpha',0.25,'edgecolor','none');


%% output
mih = [lower upper];
































