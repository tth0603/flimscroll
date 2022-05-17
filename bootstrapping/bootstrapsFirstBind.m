function mib = bootstrapsFirstBind(varargin)
% bootstrapsFirstBind(dataCells,bsNum,inarg)
% A pretty specific script. See t1p54b & t1p54e.m
%
% OUTPUT
%   a strucure with cell arrays of 
%
% USAGE
% bootstrapsFirstBind(dataCells,10,[2 0.001])
% 
% Harden 2020

datCell = varargin{1};
numOfBS = varargin{2};
inarg = varargin{3};
%here we're going to try and do everything that other bootstrap scripts do,
%but with cell arrays rather than matrices
% Initialize a shit ton of variables:
nCell = size(datCell,2); %number of data sets input
ns = zeros(nCell,1); %num of measurements for each cndn
tblCell = cell(nCell,1);
bsCell = cell(nCell,1);
intCell = cell(nCell,1);
indxCell = cell(nCell,1);
afCell = cell(nCell,1);
kCell = cell(nCell,1);
for i = 1:nCell
    intCell{i} = datCell{i}.intervals; %data cell array
    indxCell{i} = intCell{i}.*0; %make an empty cell array the size of intCell to later fill with probability indices
    ns(i) = size(intCell{i},1);
    tblCell{i} = [[1:ns(i)]' 1/ns(i)*ones(ns(i),1)]; %this is from LJF's method of picking with replacement
    bsCell{i} = zeros(ns(i),numOfBS); %attmpt to fill a cell with N x numOfBS empty mats
    afCell{i} = zeros(numOfBS,1);
    kCell{i} = zeros(numOfBS,1);
end
%here, we'll grab a BS for each cndn collectively, then fit them:

%make the BS cell - cols are BS, rose are time points
tic;
for i = 1:numOfBS
%     parfor j = 1:nCell 
    for j = 1:nCell %for each BS, cycle through each data set
        indxCell{j} = probability_steps(tblCell{j},ns(j));
        bsCell{j}(:,i) = intCell{j}(indxCell{j}); %equivalent to bs = intervals(indx); in the matrix approach
    end
end
toc
% for nBS = 10: 1.3s w/o pp; 1.3s w/pp (1.1s on 2nd iteration)
% for nBS = 1000: 112s w/o pp; 127s w/pp (123s next iteration)
%now we gotta fit them
% oriCell = datCell; %for TSing
% datCell = oriCell; % for TSing
rc = 0; %may want to change this. see t1p54b.m line 1150ish
opts1=  optimset('display','off'); %attempt to suppress max'n warnings
tic;
for i = 1:numOfBS
    for j = 1:nCell
        datCell{j}.intervals = bsCell{j}(:,i); %set the intervals for condn j to be BS i.
        outs = fminsearch('exp1_global_active_fraction_nonspecific_background_mxl',inarg,opts1,datCell(j),rc);
        afCell{j}(i) = outs(1)^2/(1 + outs(1)^2); %record these for later inspection
        kCell{j}(i) = outs(2);
    end
end
toc
%output the BS fits 
mib.afCell = afCell;
mib.kCell = kCell;
% now get the stats
afStd = zeros(nCell,1);
kStd = zeros(nCell,1);
for i = 1:nCell
    afStd(i) = std(afCell{i});
    kStd(i) = std(kCell{i});
end
mib.afStd = afStd;
mib.kStd = kStd;

























    
