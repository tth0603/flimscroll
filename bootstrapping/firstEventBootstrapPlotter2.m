function [mih,plotHandle] = firstEventBootstrapPlotter2(varargin)
%   [mih,plotHandle] = firstEventBootstrapPlotter2(bsMat,analData,plotColor,fignum)
%
% Spits out and plots the 5th & 95th percentile bootstrap vector so as to get error bars.
% pair with firstEventBootstrapGenerator2.
%
% INPUTS
%   BSmat == matstrix of bootstrapped data from firstEventBootstrapGenerator
%   analData == to get the x axis of the plot correct
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
analData = varargin{2};
colr = varargin{3};
figNum = varargin{4};

mx = ceil(max(max(bsMat)));  %mx = ceil(max(intervals)); %the longest event in data
numOfBS = size(bsMat,2);
norm = size(bsMat,1);

nukesMat = analData.nukesMat; % for the time of the sints below. from initialFractionBin2

%sort the individual bootstraps (ie the cols of bsMat)
sintsMat = zeros(ceil(max(nukesMat(:,1))) + 1,numOfBS); %plus one bc we start from zero time
for j = 1:numOfBS
    sints = [];
%     for i = 0:mx %zero bc time
    for i = 0:ceil(max(nukesMat(:,1))) % again, see initialFractionBin2.
        bsV = bsMat(:,j); %pluck out a BS
        bsV = bsV(bsV > 0); % get rid of the zeros
        logik = bsV < i;
        sints = [sints; sum(logik)];
    end
    sintsMat(:,j) = sints;
end

%normailization(?):
if norm > 150 %this is an empirical hack for the neutral condition, which has 4 reps in the resubmitted doc, but an low number of nukes per bin, so we bumped it up to 44 nukes/experiment (174 total), equivalent to the wt spacer condition. 
    norm = 176;
end
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
































    