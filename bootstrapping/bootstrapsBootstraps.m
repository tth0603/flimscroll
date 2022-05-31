function [mib,plotHandle] = bootstrapsBootstraps(varargin)
%   [mib,plotHandle] = bootstrapsBootstraps(dataMatrix,dataColumn,numOfBS,color,figNum)
% My latest and greatest bootstrap script. will spit out and plot the 5th &
% 95th percentile bootstrap vector so as to get error bars.
%
% INPUTS
%   dataMatrix == some 2D matrix with x axis values as rows
%   dataColumn == a scaler indicating the column index for your data
%   numOfBS == number of bootstrap data replicates
%
% OUTPUTS
%   mib == an N x 2 matrix of the upper and lower sints, respectively
%   plotHandle == plot handle for editing the ploted curves post hoc
%   It will also plot a shaded area replresenting the error "bars" on a
%   natural log scale
%
% Timothy Harden 2020

%initiate
dataMat = varargin{1};
col = varargin{2};
intervals = dataMat(:,col);
numOfBS = varargin{3};
n = size(dataMat,1);
mx = ceil(max(intervals)); %the longest event in data
colr = varargin{4};
figNum = varargin{5};
%make the BS mat - cols are BS, rose are time points
tbl=[[1:n]' 1/n*ones(n,1)];  %this is from LJF's method of picking with replacement
bsMat = zeros(n,numOfBS);
% tic
% parpool
parfor i = 1:numOfBS
    indx = probability_steps(tbl,n);
    bs = intervals(indx); 
    bsMat(:,i) = bs;
end
% p = gcp;
% delete(p)
% toc
% this is the slow step. timing for space:
% 1. for loop: 105s
% 2. parfor 1st time: 118s (plus par pool set up)
% 3. parfor 2nd time: 37s
% 4. parfor specifying 4 workers: 105s
%sort the individual bootstraps (ie the cols of bsMat)
sintsMat = zeros(mx + 1,numOfBS); %plus one bc we start from zero time
for j = 1:numOfBS
    sints = [];
    for i = 0:mx %zero bc time
        logik = bsMat(:,j) > i;
        sints = [sints; sum(logik)];
    end
    sintsMat(:,j) = sints;
end
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

% plot them:
figure(figNum);
% plotHandle = stairs([t t],[log(lower/n) log(upper/n)],'LineWidth',1);
% plotHandle = stairs([t t],[log(lower/n) log(upper/n)],'LineWidth',1);
% add shading
lower2 = nonzeros(lower);
upper2 = nonzeros(upper); %must get rid of zeros to make shading work bc ln leads to -inf vals
% now pad lower2 with 1's so its the same size as upper:
pad = ones(size(upper2,1) - size(lower2,1),1);
lower2 = [lower2; pad];
bnd = size(upper2,1);
lower2 = lower2';
upper2 = upper2'; %bc they need to be row v's for some reason
ts = [0:bnd - 1, fliplr(0:bnd - 1)]; %not sure why we gotta flip these things or if we really need fliplr and not just flip but whatever
shade = [log(lower2/n), fliplr(log(upper2/n))];
plotHandle = fill(ts,shade,colr); 
set(plotHandle,'facealpha',0.25,'edgealpha',0.25,'edgecolor','none');

%% output
mib = [lower upper];
































