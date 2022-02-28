function [mib,plotHandle] = firstEventBootstrap(varargin)
%   [mib,plotHandle] = firstEventBootstrap(cia,normalization,numOfBS,plotColor,fignum)
% My latest and greatest bootstrap script. will spit out and plot the 5th &
% 95th percentile bootstrap vector so as to get error bars.
%
% INPUTS
%   cia == cumulative interval array from flimscroll2
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

%initiate:
cia = varargin{1};
normalization = varargin{2};
numOfBS = varargin{3};
colr = varargin{4};
figNum = varargin{5};

%grab the first event (its length)
first = cia(cia(:,2) == -3,7); %a vector

%determine normalization
if ~isempty(normalization)
    norm = normalization;
else 
%     mx = length(first); %there is a debate in my head as to which is more
%     useful
    norm = 1;
end
% % %cumulative intervals:
% % sints = [];
% % for i = 0:ceil(max(cia(:,7)))
% %     logi = first(:) < i;
% %     sints = [sints;i sum(logi)];
% % end
% % %normailization?
% % sints(:,2) = sints(:,2)./mx;

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

% % plot TSing:
% figure(figNum);
% plotHandle = stairs([t t],[(lower') (upper')],'LineWidth',1);
% plotHandle = stairs([t t],[(lower') (upper')],'LineWidth',1);

%% output
mib = [lower upper];
































