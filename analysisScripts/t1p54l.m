% Code to reporduce all figures from Harden et al 2022 save for those
% pertaining to idle periods (i.e. Figs 3E, F, and SFig. 9; see
% t1p54r) and image acquiosition conditions (SFig 1; see t1p540)

%% load workspace
% neutral spacer eveS2
% 190712 
load t1p51o.dat -mat %this one includes so spot analysis
%reassign data:
nAnal1 = analData;
clear analData
nCia1 = nAnal1.cia;
nSpots1 = nAnal1.spotsMat; 
nNukes1 = nAnal1.nukesMat;
% 190808
load t1p58f.dat -mat
nAnal2 = analData;
clear analData
nCia2 = nAnal2.cia;
nSpots2 = nAnal2.spotsMat;
nNukes2 = nAnal2.nukesMat;
% 190809 (processes on 210806)
load t1p58k2.dat -mat
nAnal3 = analData;
clear analData
nCia3 = nAnal3.cia;
nSpots3 = nAnal3.spotsMat;
nNukes3 = nAnal3.nukesMat;
% 210910
load t1p69g2.dat -mat
nAnal4 = analData;
clear analData
nCia4 = nAnal4.cia;
nSpots4 = nAnal4.spotsMat;
nNukes4 = nAnal4.nukesMat;
% combine them
nAnal1and2 = combineAnalData(nAnal1, nAnal2);
nAnal3and4 = combineAnalData(nAnal3, nAnal4);
nAnal = combineAnalData(nAnal1and2, nAnal3and4); %another iteration bc a > 2 replicates
nCia = nAnal.cia;
nSpots = nAnal.spotsMat;
nNukes = nAnal.nukesMat;

% wt spacer
% 190726
load t1p53m.dat -mat
spaceAnal1 = analData;
clear analData
spaceCia1 = spaceAnal1.cia;
spaceSpots1 = spaceAnal1.spotsMat;
spaceNukes1 = spaceAnal1.nukesMat;
% 190815
load t1p61e2.dat -mat
spaceAnal2 = analData;
clear analData
spaceCia2 = spaceAnal2.cia;
spaceSpots2 = spaceAnal2.spotsMat;
spaceNukes2 = spaceAnal2.nukesMat;
%combined
spaceAnal = combineAnalData(spaceAnal1, spaceAnal2);
spaceCia = spaceAnal.cia;
spaceSpots = spaceAnal.spotsMat;
spaceNukes = spaceAnal.nukesMat;

% zld
% 180917
load t1p40n.dat -mat 
zldAnal1 = analData;
clear analData
zldCia1 = zldAnal1.cia;
zldSpots1 = zldAnal1.spotsMat;
zldNukes1 = zldAnal1.nukesMat;
% 191024
% load t1p59e2.dat -mat
load t1p59m.dat -mat
zldAnal2 = analData;
clear analData
zldCia2 = zldAnal2.cia;
zldSpots2 = zldAnal2.spotsMat;
zldNukes2 = zldAnal2.nukesMat;
% combine them
zldAnal = combineAnalData(zldAnal1, zldAnal2);
zldCia = zldAnal.cia;
zldSpots = zldAnal.spotsMat;
zldNukes = zldAnal.nukesMat;

% bcd
% 181010
load t1p38k.dat -mat
bcdAnal1 = analData;
clear analData;
bcdCia1 = bcdAnal1.cia;
bcdSpots1 = bcdAnal1.spotsMat;
bcdNukes1 = bcdAnal1.nukesMat;
% 191106
load t1p62e2.dat -mat
bcdAnal2 = analData;
clear analData;
bcdCia2 = bcdAnal2.cia;
bcdSpots2 = bcdAnal2.spotsMat;
bcdNukes2 = bcdAnal2.nukesMat;
% combine them
bcdAnal = combineAnalData(bcdAnal1, bcdAnal2);
bcdCia = bcdAnal.cia;
bcdSpots = bcdAnal.spotsMat;
bcdNukes = bcdAnal.nukesMat;

% dstat
% 180914
% load t1p37j.dat -mat
load t1p37m.dat -mat %same as 37j but with updated AP
dAnal1 = analData;
clear analData
dCia1 = dAnal1.cia;
dSpots1 = dAnal1.spotsMat;
dNukes1 = dAnal1.nukesMat;
% 191122
load t1p63e2.dat -mat
dAnal2 = analData;
clear analData
dCia2 = dAnal2.cia;
dSpots2 = dAnal2.spotsMat;
dNukes2 = dAnal2.nukesMat;
%combine them
dAnal = combineAnalData(dAnal1, dAnal2);
dCia = dAnal.cia;
dSpots = dAnal.spotsMat;
dNukes = dAnal.nukesMat;

% cia cell arrays for plotting from loops
ciaCell{1} = nCia; %n = 
ciaCell{2} = spaceCia; %n = 
ciaCell{3} = zldCia; %n = 
ciaCell{4} = bcdCia; %n = 
ciaCell{5} = dCia; %n = 

ciaCell1{1} = nCia1; %n = 
ciaCell1{2} = spaceCia1; %n = 
ciaCell1{3} = zldCia1; %n = 
ciaCell1{4} = bcdCia1; %n = 
ciaCell1{5} = dCia1; %n = 

ciaCell2{1} = nCia2; %n = 
ciaCell2{2} = spaceCia2; %n = 
ciaCell2{3} = zldCia2; %n = 
ciaCell2{4} = bcdCia2; %n = 
ciaCell2{5} = dCia2; %n = 

ciaCell3{1} = nCia3; % bc third replicate
ciaCell4{1} = nCia4;

% make cells limit analysis to strip center:
bins = [-0.20:0.02:0.20]; % define the bins
boiA = 5:9; % ant
boiP  = 12:16; % post (misses some events but keeps things symmetrical
allBins = 4:21; % for all events
boiC = 10:11; % for the center bins

% cia cells for plotting replicates:
% ATM we don't actually use these, but we should. It's be easier than the
% above
nCiaCell{1} = nCia1; nCiaCell{2} = nCia2; nCiaCell{3} = nCia3; nCiaCell{4} = nCia4;
spaceCiaCell{1} = spaceCia1; spaceCiaCell{2} = spaceCia2;
zldCiaCell{1} = zldCia1; zldCiaCell{2} = zldCia2;
bcdCiaCell{1} = bcdCia1; bcdCiaCell{2} = bcdCia2;
dCiaCell{1} = dCia1; dCiaCell{2} = dCia2;

% nuke Mats for plotting reps:
nNukesCell{1} = nNukes1; nNukesCell{2} = nNukes2; nNukesCell{3} = nNukes3; nNukesCell{4} = nNukes4;
spaceNukesCell{1} = spaceNukes1; spaceNukesCell{2} = spaceNukes2;
zldNukesCell{1} = zldNukes1; zldNukesCell{2} = zldNukes2;
bcdNukesCell{1} = bcdNukes1; bcdNukesCell{2} = bcdNukes2;
dNukesCell{1} = dNukes1; dNukesCell{2} = dNukes2;

% anal strucs for plotting reps:
nAnalCell{1} = nAnal1; nAnalCell{2} = nAnal2; nAnalCell{3} = nAnal3; nAnalCell{4} = nAnal4;
spaceAnalCell{1} = spaceAnal1; spaceAnalCell{2} = spaceAnal2;
zldAnalCell{1} = zldAnal1; zldAnalCell{2} = zldAnal2;
bcdAnalCell{1} = bcdAnal1; bcdAnalCell{2} = bcdAnal2;
dAnalCell{1} = dAnal1; dAnalCell{2} = dAnal2;

% center cias for fitting reps:
nCenterCia = {};
for i = 1:4
    nCenterCia{i} = binnedCiaMaker3(nAnalCell{i},bins,boiC);
end
spaceCenterCia = {};
for i = 1:2
    spaceCenterCia{i} = binnedCiaMaker3(spaceAnalCell{i},bins,boiC);
end
zldCenterCia = {};
for i = 1:2
    zldCenterCia{i} = binnedCiaMaker3(zldAnalCell{i},bins,boiC);
end
bcdCenterCia = {};
for i = 1:2
    bcdCenterCia{i} = binnedCiaMaker3(bcdAnalCell{i},bins,boiC);
end
dCenterCia = {};
for i = 1:2
    dCenterCia{i} = binnedCiaMaker3(dAnalCell{i},bins,boiC);
end

% line color set up
colr{1} = [0 0 0]; %wt blk
colr{2} = [0.6 0.6 0.6]; %wt spacer dark gray
colr{3} = [1 0.2 0.2]; %zld red
colr{4} = [0 0.447058826684952 0.74117648601532]; %bcd blue
colr{5} = [0.23 0.44 0.33]; %stat grn
colr{6} = [0.3 0.3 0.3];  %broken grey
colr{7} = [0.929 0.694 0.125]; % orange?
colr{8} = [0.494 0.184 0.556]; % purple

% legend set up 
clear tits
tits{1} = 'eve2:neutral';
tits{2} = 'eve2:wt';
tits{3} = 'eve2[Zld]:neutral';
tits{4} = 'eve2[Bcd]:neutral';
tits{5} = 'eve2[Dst]:neutral';

% figure the shortest experiment
tMxV = max(nNukes1(:,1));
tMxV = [tMxV; max(nNukes2(:,1))];
tMxV = [tMxV; max(nNukes3(:,1))];
tMxV = [tMxV; max(spaceNukes1(:,1))];
tMxV = [tMxV; max(spaceNukes2(:,1))];
tMxV = [tMxV; max(zldNukes2(:,1))];
tMxV = [tMxV; max(zldNukes2(:,1))];
tMxV = [tMxV; max(bcdNukes2(:,1))];
tMxV = [tMxV; max(bcdNukes2(:,1))];
tMxV = [tMxV; max(dNukes2(:,1))];
tMxV = [tMxV; max(dNukes2(:,1))];
tMx = min(tMxV);

analCell1{1} = nAnal1;
analCell2{1} = nAnal2;
analCell3{1} = nAnal3; % bc thrid rep
analCell4{1} = nAnal4;
analCell1{2} = spaceAnal1;
analCell2{2} = spaceAnal2;
analCell1{3} = zldAnal1;
analCell2{3} = zldAnal2;
analCell1{4} = bcdAnal1;
analCell2{4} = bcdAnal2;
analCell1{5} = dAnal1;
analCell2{5} = dAnal2;

analCell{1} = nAnal;
analCell{2} = spaceAnal;
analCell{3} = zldAnal;
analCell{4} = bcdAnal;
analCell{5} = dAnal;
% make the center nuke cell arrays:
centerCia = cell(1,5);centerCia1 = cell(1,5);centerCia2 = cell(1,5);centerCia3 = cell(1,5);centerCia4 = cell(1,5);
%take care of neutral first bc of the extra reps
centerCia{1} = binnedCiaMaker3(analCell1{1},analCell2{1},analCell3{1},analCell4{1},bins,boiC);
centerCia1{1} = (binnedCiaMaker3(analCell1{1},bins,boiC));
centerCia2{1} = (binnedCiaMaker3(analCell2{1},bins,boiC));
centerCia3{1} = (binnedCiaMaker3(analCell3{1},bins,boiC));
centerCia4{1} = (binnedCiaMaker3(analCell4{1},bins,boiC));
% now the others:
for i = 2:5
    centerCia{i} = (binnedCiaMaker3(analCell1{i},analCell2{i},bins,boiC));
    centerCia1{i} = (binnedCiaMaker3(analCell1{i},bins,boiC));
    centerCia2{i} = (binnedCiaMaker3(analCell2{i},bins,boiC));
end
%% nuke binning stuff
% first, figure out how many nukes are in each bin, as this can vary from
% expreiment to expreiment, but in reality should be more or less constant:
% here is a matrix of the number of nuclei in each bin:
totalBinMat = [];
for i = 1:5
    out = positionDistribution(analCell1{i},bins);
    totalBinMat = [totalBinMat out.totalBin];
    out = positionDistribution(analCell2{i},bins);
    totalBinMat = [totalBinMat out.totalBin];
    try
        out = positionDistribution(analCell3{i},bins);
        totalBinMat = [totalBinMat out.totalBin];
        out = positionDistribution(analCell4{i},bins);
        totalBinMat = [totalBinMat out.totalBin];
    end
end
% its not fair to count the zeros in determining an average, so convert
% those to NANs
totalBinMat(totalBinMat == 0) = NaN;
totalBinMean = nanmean(totalBinMat,2); % we could sub this in for out.totalBin values below. It doesn't change things much and is maybe more honest? --211001
medNukeNum = median(nansum(totalBinMat)); %maybe we use this for a normalized spots v time plot?
medCenterNukeNum = median(nansum(totalBinMat(10:11,:))); % use this for spots v time?
% now get a matrix of all the bins with nukes:
eventBinMat = [];
for i = 1:5
    out = positionDistribution(analCell1{i},bins);
    eventBinMat = [eventBinMat out.eventBin];
    out = positionDistribution(analCell2{i},bins);
    eventBinMat = [eventBinMat out.eventBin];
    try
        out = positionDistribution(analCell3{i},bins);
        eventBinMat = [eventBinMat out.eventBin];
        out = positionDistribution(analCell4{i},bins);
        eventBinMat = [eventBinMat out.eventBin];
    end
end
sum(eventBinMat);

%% example trace (from pre-review version)
figNum = 30;
nukeVect = 348;
specifiedFullTraceArrayFli(spaceAnal1,nukeVect,figNum);
%% a number of traces (for resubmission, SFig 2)
% wt:
% to select specific nukes:
fignum = 100;
wtv = [228 236 196 177 162 168 205 207 90];
% traceArrayFli(nNukes1,200) % to look at all the traces
%%randomly pick traces:
% nukesMat = nNukes1;
% nukeV = unique(nukesMat(:,3)); %all nuke numbers
% NOI = [];
% newNukesMat = [];
% for i = nukeV'
%     logi = nukesMat(:,3) == i;
%     mat = nukesMat(logi,:);
%     if sum(mat(:,4)) > 0
%         NOI = [NOI; i];
%         newNukesMat = [newNukesMat; mat];  %if the current nuke shows even one spot, we add that data to a new mat
%     end
% end
% newNukeV = unique(newNukesMat(:,3));
% n = length(newNukeV);
% tbl=[[1:n]' 1/n*ones(n,1)];  %this is from LJF's method of picking with replacement
% indx = probability_steps(tbl,n); 
% wtv = newNukeV(indx(1:9)'); 
specifiedFullTraceArrayFli(nAnal1,wtv,fignum);
%%test out getting rid of one frame dropouts:
%%nvm. I don't think this is an issue. Or it will just create another issue. as is, the way I plot the traces already doesn't display single frame drop outs. 
%% save:
fp = 'C:\Users\thard\Dropbox\hardenTFfunction2020\figures\traces\neutral190712\';
for i = 1:9
    fn = sprintf('%snuke%d',fp,wtv(i));
    f = fignum -1 + i; % -1 bc indexing error
    saveas(figure(f),fn,'svg');
end

%% space traces - SFig 2
fignum = 110;
spacev = [315 317 319 321 290 306 283 271 143]; % 
% spacev = [166 222 226 283 284 290 315 319 348];
% traceArrayFli(spaceNukes,150);
% specifiedTraceArrayFli(spaceNukes,spacev,100);
%%randomly pick traces:
% nukesMat = spaceNukes1;
% nukeV = unique(nukesMat(:,3)); %all nuke numbers
% NOI = [];
% newNukesMat = [];
% for i = nukeV'
%     logi = nukesMat(:,3) == i;
%     mat = nukesMat(logi,:);
%     if sum(mat(:,4)) > 0
%         NOI = [NOI; i];
%         newNukesMat = [newNukesMat; mat];  %if the current nuke shows even one spot, we add that data to a new mat
%     end
% end
% newNukeV = unique(newNukesMat(:,3));
% n = length(newNukeV);
% tbl=[[1:n]' 1/n*ones(n,1)];  %this is from LJF's method of picking with replacement
% indx = probability_steps(tbl,n); 
% spacev = newNukeV(indx(1:9)'); 
specifiedFullTraceArrayFli(spaceAnal1,spacev,fignum);
%% save:
fp = 'C:\Users\thard\Dropbox\hardenTFfunction2020\figures\traces\wtSpacer190726\';
for i = 1:9
    fn = sprintf('%snuke%d',fp,spacev(i));
    f = fignum -1 + i; % -1 bc indexing error
    saveas(figure(f),fn,'svg');
end

%% spots v T - Fig 1D
%Importantly, here we normalize by the median number of nukes in all
%experiments. Run the nukes binning section for that
% neutral:
figure(35);
clear mat mat2 mn mx
%get all rose with a spot:
mat = nNukes(nNukes(:,4) > 0,:); %col 4 is binary
size(unique(mat(:,3)))
%get frame range:
mx = max(mat(:,2)); %frame indexing is more useful that time
svt = []; %spots v time
for i = 1:mx
    mat2 = mat(mat(:,2) == i,:); %find all rose for this frame (logical indexing)
    if ~isempty(mat2) %if we have nukes with spots, mark the time and sum up the num of spots
        svt = [svt;mat2(1,1) size(mat2,1)]; %[time numOfSpots]
    end
end
hold on;plot(svt(:,1),svt(:,2)./(2*medNukeNum),'Color',colr{1},'LineStyle','-','Marker','.','MarkerSize',10);shg % x 2 bc double the reps

% space:
clear mat mat2 mn mx
%get all rose with a spot:
mat = spaceNukes(spaceNukes(:,4) > 0,:); %col 4 is binary
%get frame range:
mx = max(mat(:,2)); %frame indexing is more useful that time
svt = []; %spots v time
for i = 1:mx
    mat2 = mat(mat(:,2) == i,:); %find all rose for this frame (logical indexing)
    if ~isempty(mat2) %if we have nukes with spots, mark the time and sum up the num of spots
        svt = [svt;mat2(1,1) size(mat2,1)]; %[time numOfSpots]
    end
end
hold on;plot(svt(:,1),svt(:,2)./medNukeNum,'Color',colr{2},'LineStyle','-','Marker','.','MarkerSize',10);shg

%  zld
clear mat mat2 mn mx
%get all rose with a spot:
mat = zldNukes(zldNukes(:,4) > 0,:); %col 4 is binary
%get frame range:
mx = max(mat(:,2)); %frame indexing is more useful that time
svt = []; %spots v time
for i = 1:mx
    mat2 = mat(mat(:,2) == i,:); %find all rose for this frame (logical indexing)
    if ~isempty(mat2) %if we have nukes with spots, mark the time and sum up the num of spots
        svt = [svt;mat2(1,1) size(mat2,1)]; %[time numOfSpots]
    end
end
hold on;plot(svt(:,1),svt(:,2)./medNukeNum,'Color',colr{3},'LineStyle','-','Marker','.','MarkerSize',10);shg

%  bcd
clear mat mat2 mn mx
%get all rose with a spot:
mat = bcdNukes(bcdNukes(:,4) > 0,:); %col 4 is binary
%get frame range:
mx = max(mat(:,2)); %frame indexing is more useful that time
svt = []; %spots v time
for i = 1:mx
    mat2 = mat(mat(:,2) == i,:); %find all rose for this frame (logical indexing)
    if ~isempty(mat2) %if we have nukes with spots, mark the time and sum up the num of spots
        svt = [svt;mat2(1,1) size(mat2,1)]; %[time numOfSpots]
    end
end
hold on;plot(svt(:,1),svt(:,2)./medNukeNum,'Color',colr{4},'LineStyle','-','Marker','.','MarkerSize',10);shg

% stat
clear mat mat2 mn mx
%get all rose with a spot:
mat = dNukes(dNukes(:,4) > 0,:); %col 4 is binary
%get frame range:
mx = max(mat(:,2)); %frame indexing is more useful that time
svt = []; %spots v time
for i = 1:mx
    mat2 = mat(mat(:,2) == i,:); %find all rose for this frame (logical indexing)
    if ~isempty(mat2) %if we have nukes with spots, mark the time and sum up the num of spots
        svt = [svt;mat2(1,1) size(mat2,1)]; %[time numOfSpots]
    end
end
hold on;plot(svt(:,1),svt(:,2)./medNukeNum,'Color',colr{5},'LineStyle','-','Marker','.','MarkerSize',10);shg
set(gca,'Box',true,'FontSize',16);
set(gca,'YLim',[0 0.65],'XLim',[0 4400]);
% title('Total spots v. time');legend(tits);
% 210930: the problem with spots v t is that is no longer a fair
% comparison, now that we have 4 embryos for neutral

%% spots v time in the center of the stripe: (211005 - unused)
%Importantly, here we normalize by the median number of nukes in all
%experiments. Run the nukes binning section for that
% neutral:
figure(4);
clear mat mat2 mn mx
for i = 1:5
    bin1 = 10; bin2 = 11; % these are the two indicies within BINS that constitute the center of the stripe
    % for rep 1:
    % get analdata & nuke data:
    anal = analCell1{i};
    nukesMat = anal.nukesMat;
    % get nukes with events:
    nukeEvent = nukesMat(nukesMat(:,4) > 0,:); 
    %get frame range:
    mx = max(nukeEvent(:,2));
    % bin stuff:    
    out = positionDistribution(anal,bins);
    nukePos = out.eventPos; % the nuke nums and their binning (col 3)
    % now get the center nukes:
    centerNukes = nukePos(nukePos(:,3) == bin1 | nukePos(:,3) == bin2,:);
    % now take only center nukes from the nukesMat:
    nuke = [];
    for j = 1:length(centerNukes)
        mat = nukeEvent(nukeEvent(:,3) == centerNukes(j),:);
        nuke = [nuke;mat];
    end
    svt1 = [];
    for j = 1:mx % iterate on frames
        frameNukes = nuke(nuke(:,2) == j,:);
        if ~isempty(frameNukes) % if we have nukes here, we check to see if theyre in the center
            svt1 = [svt1;frameNukes(1,2) size(frameNukes,1)];
        end
    end
% for second rep
    % get analdata & nuke data:
    anal = analCell2{i};
    nukesMat = anal.nukesMat;
    % get nukes with events:
    nukeEvent = nukesMat(nukesMat(:,4) > 0,:); 
    %get frame range:
    mx = max(nukeEvent(:,2));
    % bin stuff:
    out = positionDistribution(anal,bins);
    nukePos = out.eventPos; % the nuke nums and their binning (col 3)
    % now get the center nukes:
    centerNukes = nukePos(nukePos(:,3) == bin1 | nukePos(:,3) == bin2,:);
    % now take only center nukes from the nukesMat:
    nuke = [];
    for j = 1:length(centerNukes)
        mat = nukeEvent(nukeEvent(:,3) == centerNukes(j),:);
        nuke = [nuke;mat];
    end
    svt2 = [];
    for j = 1:mx % iterate on frames
        frameNukes = nuke(nuke(:,2) == j,:);
        if ~isempty(frameNukes) % if we have nukes here, we check to see if theyre in the center
            svt2 = [svt2;frameNukes(1,2) size(frameNukes,1)];
        end
    end
    % combine the svts:
    mxFrames = max([svt1(:,1);svt2(:,1)]);
    svt = [];
    for j = 1:mxFrames
        if any(svt1(:,1) == j) && any(svt2(:,1) == j) % if both svts have spots in frame i
            n1 = svt1(svt1(:,1) == j,2);
            n2 = svt2(svt2(:,1) == j,2);
            svt = [svt;j n1 + n2]; %add that number to the number of spots in svt1
        elseif any(svt1(:,1) == j) % if only svt1 has spots in this frame
            svt = [svt;svt1(svt1(:,1) == j,:)];
        elseif any(svt2(:,1) == j) % if only svt2 has spots in this frame
            svt = [svt;svt2(svt2(:,1) == j,:)];
        end
    end
    
hold on;plot(svt(:,1),svt(:,2),'Color',colr{i},'LineStyle','-','Marker','.','MarkerSize',10);shg    
end
set(gca,'Box',true,'FontSize',16);
ylabel('number of spots');xlabel('time (s)');title('Total spots v. time');legend(tits);

% for some reason the above is not picking out center spots. Or it gives
% the same result as the spots v t above. I am not sure why this is. maybe
% try limiting the nukes to center nukes at the start of the routine. 

%also, if you ever do get this working, you have to convert frames back to
%time at the end

%% First passage: global gamma fit: global rate, active fraction. fig. 2B. 
% center nukes, norm'd to one, Fit the rate (scale
% parameter,B) constant globally:

% plot data first:
fignum = 17;
NtV = zeros(1,5);
% take care of neutral 1st (with its extra reps):
[out,p] = initialFractionBin2(analCell1{1},analCell2{1},analCell3{1},analCell4{1},bins,boiC,[],30,3000,fignum);
% see initialFractionBin3 for exceptions to normalization that come from
% analyzing the two new reps of neutral
% [out,p] = initialFractionBin2(analCell1{1},analCell2{1},bins,boiC,[],30,3000,fignum); % the 1st submission method, with only two neutral reps
set(p,'Color',colr{1},'LineWidth',1);
hold on
NtV(1) = out.totNukes;
% now the rest:
for i = 2:5
    [out,p] = initialFractionBin2(analCell1{i},analCell2{i},bins,boiC,[],30,3000,fignum);
    set(p,'Color',colr{i},'LineWidth',1);
    hold on
    NtV(i) = out.totNukes;
end
set(gca,'YLim',[0 1.01],'Xlim',[0 4000],'Box',true,'FontSize',16); 
ylabel('Fraction of nuclei')
xlabel('Time (s)')
% fit:
inarg = [200 12 4 4 8 8 0.4    0.9    0.9    0.7    0.9]; 
gFit4 = fminsearch('globalGammaFit4',inarg,[],centerCia,NtV);
% add fits to plot:
xs =0:4000;
rate = gFit4(1);
stepsV = gFit4(2:6);
AfV = gFit4(7:end);
for i = 1:5
    gCDF = AfV(i)*gamcdf(xs,stepsV(i),rate);
    hold on;plot(xs,gCDF,'--','Color',colr{i},'LineWidth',0.75);shg
end
% here are the fit parameters saved as the variable gFit4:
% save C:\Users\thard\Dropbox\matlab\data\t1p54k4.dat gFit4

%% save
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(gca,[fp 'firstPassageCurve&Fits220328'],'svg'); %this is for the resubmission. Identical to the original, but with the two extra neutral reps included. 
% the normalization for the neutral conditions are hacky, but well
% justified (see initialFractionBin2 ~line 112; firstEventBootstrapPlotter2
% ~line 48). 10000 bootstraps. 

%% first bind BS - for Fig. 2B
% generate bootstraps
bs = 10000; %211005 - i used 10000 BSs (~4.8 hours on my laptop) for the re-submission. excessive.
bsCell = {};
% neutral first:
tic;
bsCell{1} = firstEventBootstrapGenerator2(analCell1{1},analCell2{1},analCell3{1},analCell4{1},bins,boiC,bs);
toc
% now the rest:
for i = 2:5
    tic;
    bsCell{i} = firstEventBootstrapGenerator2(analCell1{i},analCell2{i},bins,boiC,bs);
    toc
end

% first bind fit error analysis
inarg = [200 12 4 4 8 8 0.4    0.9    0.9    0.7    0.9]; 
tic;
fitErrStruc = firstEventBootstrapFitter2(bsCell,inarg);
toc
% saved so you dont ahve to run again:
% save C:\Users\thard\Dropbox\matlab\data\t1p54k1.dat fitErrStruc bsCell
% 211006 here are 10000 BSs for the resubmission figure:
% save C:\Users\thard\Dropbox\matlab\data\t1p54k2.dat fitErrStruc bsCell

%% 211008 enhanced error analysis - for FIg. 2B
% bc the 10000 bootstraps look fin on the plot, but have have some insane
% outliers in the fitting, let's redo the error using only the first 97%
% of bootstraps. The outliers do not represent sensical fit
% parameters:
% bs = 10000;
fitMat = fitErrStruc.fitMat;
fitMatSort = sort(fitMat);
% pluck out the 5th & 95th indicies:
upperLimit = ceil(0.97*bs);
fitMat90 = fitMatSort(1:upperLimit,:);
fitErrStruc.mean = mean(fitMat90);
fitErrStruc.std = std(fitMat90);
% 211008 are the bootstraps and the revised error analysis (these are the
% numbers that made it into the v3 version of the resubmission tables):
% save C:\Users\thard\Dropbox\matlab\data\t1p54k3.dat fitErrStruc bsCell


%% plot bootstraps - Fig. 2B
fignum = 22; 
% load C:\Users\thard\Dropbox\matlab\data\t1p54k1.dat -mat % original submision
% load C:\Users\thard\Dropbox\matlab\data\t1p54k2.dat -mat % re-submission 211006  (this takes a long time to plot -- 7 minutes on my laptop)
tic;
for i = 1:5
    [~,p] = firstEventBootstrapPlotter2(bsCell{i},analCell{i},colr{i},fignum);
    hold on;
end
set(gca,'YLim',[0 1.01],'Xlim',[0 4000],'Box',true,'FontSize',16); 
ylabel('Fraction of nuclei')
xlabel('Time (s)')
toc;
%% save
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(gca,[fp 'firstPassageBS220328'],'svg'); % pair with firstPassageCurve&Fits211006

%% surv Freq of active trxn, center (unused)
% -- choose this one
timestep = 15;
fignum = 30;
activationParamsCenterDropout = zeros(5,3);
freqCenterDropout = zeros(5,1);
eventNumberCenterDropout = zeros(5,1);
lowTimeCenterDropout = zeros(5,1);
expLengthMat = zeros(5,4);
meanExpLength = 3.2730e+03; % I got this from averaging all ten exp lengths computed in frequency_dwellFli2
for i = 1:5
%     out = frequency_dwellFli3(dropout(centerCia1{i}),dropout(centerCia2{i}),timestep,colr{i},fignum,[0.9 0.01 0.001],[]); % 
%     expLengthMat(i,1:2) = out.expLength;
    try
%         out = frequency_dwellFli5(dropout(centerCia1{i}),dropout(centerCia2{i}),dropout(centerCia3{i}),dropout(centerCia4{i}),timestep,colr{i},fignum,[0.8 0.1 0.01],[]); % dropout, center
        out = frequency_dwellFli5((centerCia1{i}),(centerCia2{i}),(centerCia3{i}),(centerCia4{i}),timestep,colr{i},fignum,[0.8 0.1 0.01],[]); % center
%         out = frequency_dwellFli5(dropout(ciaCell1{i}),dropout(ciaCell2{i}),dropout(ciaCell3{i}),dropout(ciaCell4{i}),timestep,colr{i},fignum,[0.8 0.1 0.01],[]); % dropout, whole stripe
            expLengthMat(i,:) = out.expLength;
    catch
%         out = frequency_dwellFli5(dropout(centerCia1{i}),dropout(centerCia2{i}),timestep,colr{i},fignum,[0.8 0.1 0.01],[]); % dropout, center
        out = frequency_dwellFli5((centerCia1{i}),(centerCia2{i}),timestep,colr{i},fignum,[0.8 0.1 0.01],[]); % center
        expLengthMat(i,1:2) = out.expLength;
    end
    activationParamsCenterDropout(i,:) = [out.a out.k1 out.k2];
    freqCenterDropout(i) = 1/out.totalFrequency;
    eventNumberCenterDropout(i) = out.eventNumber;
    lowTimeCenterDropout(i) = out.lowTimes;
end
set(gca,'YLim',[1e-5 1.3e-2],'XLim',[0 2300]);

%% active trxn survival & fit; Fig. 3B
fignum = 31;
numOfBS = 10000;
flag = 1;

inarg = [0.8 0.1 0.01];
activationParamsCenter = zeros(5,3);
for i = 1:5
    [out,p] = survivalPlotFromVector((centerCia{i}),6,2,1,3600,fignum,inarg,colr{i}); %
    set(p,'Color',colr{i});
    if flag == 1 
        [~,pbs] = bootstrapsBootstraps((centerCia{i}),6,numOfBS,colr{i},fignum);
    end
    hold on;
    activationParamsCenter(i,:) = [out.a out.k1 out.k2];
end
set(gca,'Box',true,'FontSize',16,'YLim',[-6.5 0],'XLim',[0 2300]);

%% save this boot strap plot:
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
% saveas(gca,[fp 'survFreqCurveFit&err220328'],'svg'); % this includes the new, full analysis of the zld rep from 191024
% saveas(gca,[fp 'activeTrxnSurvFit&err220402'],'svg'); % this includes dropouts
saveas(gca,[fp 'activeTrxnSurvFit&err220412'],'svg'); % this is the survival version of the plot

%% surv Freq, center, bootstraps (unused)
timestep = 15;
fignum = 33;
constructs = 1:5;
tic;
for i = constructs
%     [~,pbs] = bootstrapsFrequency(dropout(centerCia{i}),6,1000,freqCenterDropout(i),colr{i},fignum);
    [~,pbs] = bootstrapsFrequency(centerCia{i},6,10000,freqCenterDropout(i),colr{i},fignum);
    hold on
end
toc
% this will take about 5 min to run on your lab machine with bs = 10000
%
for i = constructs
    try
%         out = frequency_dwellFli5(dropout(centerCia1{i}),dropout(centerCia2{i}),dropout(centerCia3{i}),dropout(centerCia4{i}),timestep,colr{i},fignum,[0.9 0.01 0.001],[]); % 
        out = frequency_dwellFli5((centerCia1{i}),(centerCia2{i}),(centerCia3{i}),(centerCia4{i}),timestep,colr{i},fignum,[0.9 0.01 0.001],[]); % 
    catch
%         out = frequency_dwellFli5(dropout(centerCia1{i}),dropout(centerCia2{i}),timestep,colr{i},fignum,[0.9 0.01 0.001],[]); % 
        out = frequency_dwellFli5((centerCia1{i}),(centerCia2{i}),timestep,colr{i},fignum,[0.9 0.01 0.001],[]); % 
    end
    hold on
end
set(gca,'YLim',[1e-5 1.3e-2],'XLim',[0 2300]);
%% save this boot strap plot:
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
% saveas(gca,[fp 'survFreqCurveFit&err220328'],'svg'); % this includes the new, full analysis of the zld rep from 191024
saveas(gca,[fp 'survFreqCurveFit&err220401'],'svg'); % this includes dropouts

%% active param errors
inarg = [1 3600 0.9 0.01 0.001]; %[tm tx a k1 k2]
% parpool % use this to activate interactive session
% delete(gcp('nocreate')) %use this to end interactive session
tic;
constructs = 1:5;
activeCenterFits = {};
activationParamsCenterDropoutStd1000 = zeros(5,3);
% parfor i = constructs
for i = constructs
    bsFits = bootstrapsBootstrapsFit((centerCia{i}),6,numOfBS,2,inarg);
    activeCenterFits{i} = bsFits;
    activationParamsCenterDropoutStd1000(i,:) = bsFits.std;
end
toc
% should take about 8 min on your lab machine with bs = 10000

%% save the BSs:
save C:\Users\tth12\Dropbox\matlab\data\activeParamErra220401.dat activeCenterFits activationParamsCenterDropoutStd1000

% NB: I've looked at fitting both with rates and taus. They are for the
% most part the same, with the rate fitting a bit better behaved (smaller,
% more reasonable error bars). The one difference in neutral, that gets a
% bit crazy with the taus fit. 201201

%% table. active params. center, TABLE 2
% set up table
Reporter = tits';
% rates:
% allocate:
A = zeros(5,1);k1 = A; k2 = A; %allocate space, amke things col vectors
Aerr = A; k1err = A; k2err = A;
errMat = activationParamsCenterDropoutStd1000;
for i = 1:5
    A(i) = round(activationParamsCenter(i,1),2);
    Aerr(i) = round(errMat(i,1),2);
    k1(i) = (activationParamsCenter(i,2));
    k1err(i) = round(errMat(i,2),3);
    k2(i) = (activationParamsCenter(i,3));
    k2err(i) = round(errMat(i,3),4);
end
% onTimes = table(Reporter,A,Aerr,k1,k1err,k2,k2err)
% taus:
t1 = 1./k1;
t2 = 1./k2;
%propagate errror:
t1err = round(1./k1.^2.*k1err);
t2err = round(1./k2.^2.*k2err,-1);
% % round
t1 = round(t1);
t2 = round(t2,-1);
activeParamsCenterDropout = table(Reporter,A,Aerr,t1,t1err,t2,t2err)

% NB: 211014 I updated v3 of the tables (in C:\Users\tth12\Dropbox\hardenTFfunction2020\210226cellSystems) with these values 

%% bar chart. active params. center, Fig 3B
ax201 = figure(201);
bar(1,A(1),0.4,'FaceColor',colr{1},'EdgeColor','none'); %initiate the bar chart first
for i = 2:5
    hold on; bar(i,A(i),0.4,'FaceColor',colr{i},'EdgeColor','none');
end
xticks([1 2 3 4 5]);
set(gca,'xticklabels',tits(1:5),'yLim',[0 1],'Box',true,'FontSize',16,'xLim',[0.5 5.5]);
xticklabel_rotate([],45,[],'FontSize',16);
hold on;
er = errorbar(1,A(1),Aerr(1),Aerr(1));
er.Color = colr{2};    %make the black neutral grey                       
er.LineStyle = 'none'; 
for i = 2:5
    hold on;
    er = errorbar(i,A(i),Aerr(i),Aerr(i));
    er.Color = colr{1};    % make all err bars black                       
    er.LineStyle = 'none'; 
end
set(gca,'Box',true);
% k1's -- full distribution
ax202 = figure(202);
bar(1,t1(1),0.4,'FaceColor',colr{1},'EdgeColor','none');
for i = 2:5
    hold on; bar(i,t1(i),0.4,'FaceColor',colr{i},'EdgeColor','none');
end
xticks([1 2 3 4 5]);
set(gca,'xticklabels',tits(1:5),'yLim',[0 70],'Box',true,'FontSize',16,'xLim',[0.5 5.5]);
xticklabel_rotate([],45,[],'FontSize',16);
hold on;
er = errorbar(1,t1(1),t1err(1),t1err(1));
er.Color = colr{2};                           
er.LineStyle = 'none'; 
for i = 2:5
    hold on;
    er = errorbar(i,t1(i),t1err(i),t1err(i));
    er.Color = colr{1};                           
    er.LineStyle = 'none'; 
end
set(gca,'Box',true);
% k2's -- full distribution
ax203 = figure(203);
bar(1,t2(1),0.4,'FaceColor',colr{1},'EdgeColor','none');
for i = 2:5
    hold on; bar(i,t2(i),0.4,'FaceColor',colr{i},'EdgeColor','none');
end
xticks([1 2 3 4 5]);
set(gca,'xticklabels',tits(1:5),'yLim',[0 430],'Box',true,'FontSize',16,'xLim',[0.5 5.5]);
xticklabel_rotate([],45,[],'FontSize',16);
hold on;
er = errorbar(1,t2(1),t2err(1),t2err(1));
er.Color = colr{2};                           
er.LineStyle = 'none'; 
for i = 2:5
    hold on;
    er = errorbar(i,t2(i),t2err(i),t2err(i));
    er.Color = colr{1};                           
    er.LineStyle = 'none'; 
end
set(gca,'Box',true);
%% save
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(ax201,[fp 'freqBarAs220401'],'svg');
saveas(ax202,[fp 'freqBarTau1220401'],'svg');
saveas(ax203,[fp 'freqBarTau220401'],'svg');

% NB: updated and saved 211014

%% active trxn surv , EARLY, center, Fig. 3C
%used in paper figs
fignum = 35;
flag = 1; % 1 for bootstraps, 0 for not
inarg = [0.9 0.01 0.001];
numOfBS = 10000;
split = 1/5; %Sets the divide btwn the two bins set to 1/2 for 50/50, 1/3 for first third in "early" bin, etc
timestep = 15;
onHalf1cell1 = {};onHalf1cell2 = {};freqCenterDropoutEarly = [];
activationParamsCenterDropoutEarly = zeros(5,3);
freqCenterDropoutEarly = zeros(5,1);
eventNumberCenterDropoutEarly = zeros(5,1);
lowTimeCenterDropoutEarly = zeros(5,1);expLengthMatEarly = [];
constructs = 1:5;
for i = constructs
    % treat each rep indie:
    mat1 = centerCia1{i};
    L = size(mat1,1);
    ts = sort(mat1(:,7));
    t50 = ts(floor(L*split));
    logi = mat1(:,7) <= t50;
    half1mat = mat1(logi,:);
%     onHalf1cell1{i} = dropout(half1mat);
    onHalf1cell1{i} = (half1mat);
    num = sum(logi);
    denom = size(half1mat);
    % the other rep:
    mat2 = centerCia2{i};
    L = size(mat2,1);
    ts = sort(mat2(:,7));
    t50 = ts(floor(L*split));
    logi = mat2(:,7) <= t50;
    half1mat = mat2(logi,:);
%     onHalf1cell2{i} = dropout(half1mat);
    onHalf1cell2{i} = (half1mat);
    % maybe there are other reps:
    try
        % rep3
        mat3 = centerCia3{i};
        L = size(mat3,1);
        ts = sort(mat3(:,7));
        t50 = ts(floor(L*split));
        logi = mat3(:,7) <= t50;
        half1mat = mat3(logi,:);
%         onHalf1cell3{i} = dropout(half1mat);
        onHalf1cell3{i} = (half1mat);
        % rep4:
        mat4 = centerCia4{i};
        L = size(mat4,1);
        ts = sort(mat4(:,7));
        t50 = ts(floor(L*split));
        logi = mat4(:,7) <= t50;
        half1mat = mat4(logi,:);
%         onHalf1cell4{i} = dropout(half1mat);
        onHalf1cell4{i} = (half1mat);
        % for ferquency:
% %         [out,p] = frequency_dwellFli2(onHalf1cell1{i},onHalf1cell2{i},onHalf1cell3{i},onHalf1cell4{i},timestep,colr{i},fignum,inarg,[]);
%         out = frequency_dwellFli5(onHalf1cell1{i},onHalf1cell2{i},onHalf1cell3{i},onHalf1cell4{i},timestep,colr{i},fignum,inarg,[]); % use this for freq
% %         expLengthMatEarly(i,:) = out.expLength;
        % for survival
        [out,p] = survivalPlotFromVector([onHalf1cell1{i};onHalf1cell2{i};onHalf1cell3{i};onHalf1cell4{i}],6,2,1,3600,fignum,inarg,colr{i}); %
        set(p,'Color',colr{i});
    catch
        % for ferquency:
% %         [out,p] = frequency_dwellFli2(onHalf1cell1{i},onHalf1cell2{i},timestep,colr{i},fignum,inarg,[]);
%         out = frequency_dwellFli5(onHalf1cell1{i},onHalf1cell2{i},timestep,colr{i},fignum,inarg,[]); % use this for freq
% %         expLengthMatEarly(i,:) = out.expLength; 
        % for survival
        [out,p] = survivalPlotFromVector([onHalf1cell1{i};onHalf1cell2{i}],6,2,1,3600,fignum,inarg,colr{i}); %
        set(p,'Color',colr{i});
    end
    activationParamsCenterDropoutEarly(i,:) = [out.a out.k1 out.k2];
    % for frequency:
%     freqCenterDropoutEarly(i) = 1/out.totalFrequency;
%     eventNumberCenterDropoutEarly(i) = out.eventNumber;
%     lowTimeCenterDropoutEarly(i) = out.lowTimes;
    hold on
    if flag == 1
        try
            % frequency:
%             [~,pbs] = bootstrapsFrequency([onHalf1cell1{i};onHalf1cell2{i};onHalf1cell3{i};onHalf1cell4{i}],6,numOfBS,freqCenterDropoutEarly(i),colr{i},fignum);
            % survival:
            [~,pbs] = bootstrapsBootstraps([onHalf1cell1{i};onHalf1cell2{i};onHalf1cell3{i};onHalf1cell4{i}],6,numOfBS,colr{i},fignum);
        catch
            % fdrequency:
%             [~,pbs] = bootstrapsFrequency([onHalf1cell1{i};onHalf1cell2{i}],6,numOfBS,freqCenterDropoutEarly(i),colr{i},fignum);
            % survial:
            [~,pbs] = bootstrapsBootstraps([onHalf1cell1{i};onHalf1cell2{i}],6,numOfBS,colr{i},fignum);
        end
        hold on
    end
end
% For freque3ncy:
% set(gca,'YLim',[2e-4 0.05],'XLim',[0 2300]); % all cndns
% for survival:
set(gca,'Box',true,'FontSize',16,'YLim',[-4.8 0],'XLim',[0 2300]);

%% save
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
% saveas(gca,[fp 'survFreqCurveFit&errEarly211014'],'svg'); %saved with 10000 BSs
% saveas(gca,[fp 'survFreqCurveFit&errEarly220401'],'svg'); %this is with dropouts
saveas(gca,[fp 'survFreqCurveFit&errEarly220412'],'svg'); % the survival version

%% error surv Freq, EARLY, center, Fig. 3C
fitParams = [1 3600 0.9 0.01 0.001]; %[tm tx a k1 k2]
bsNum = 10000;
% parpool
tic;
activationParamsCenterDropoutEarlyStd1000 = [];activeCenterEarlyFits = {};activeCenterEarlyFitsTaus = {};
activationParamsCenterDropoutEarlyStdTaus1000 = [];
constructs = 1:5;
for i = 1:5
    try
        mat = [onHalf1cell1{i}; onHalf1cell2{i}; onHalf1cell3{i}; onHalf1cell4{i}];
    catch
        mat = [onHalf1cell1{i}; onHalf1cell2{i}];
    end
%     activeCenterEarlyFitsTaus{i} = bootstrapsBootstrapsFitTaus(mat,6,bsNum,2,fitParams);
%     activationParamsCenterDropoutEarlyStdTaus1000(i,:) = activeCenterEarlyFitsTaus{i}.std;
    activeCenterEarlyFits{i} = bootstrapsBootstrapsFit(mat,6,bsNum,2,fitParams);
    activationParamsCenterDropoutEarlyStd1000(i,:) = activeCenterEarlyFits{i}.std;
end
toc

%% save the Bootstraps:
save C:\Users\thard\Dropbox\matlab\data\freqSurvEarlyErrs211015.dat activeCenterEarlyFitsTaus activationParamsCenterDropoutEarlyStdTaus1000 activeCenterEarlyFits activationParamsCenterDropoutEarlyStd1000

%% table surv Freq, EARLY, center, TABLE 2
Reporter = tits';
% rates:
% allocate:
EA = zeros(5,1);Ek1 = EA; Ek2 = EA; %allocate space, amke things col vectors
EAerr = EA; Ek1err = EA; Ek2err = EA;
errMat = activationParamsCenterDropoutEarlyStd1000;
for i = 1:5
    EA(i) = round(activationParamsCenterDropoutEarly(i,1),2);
    EAerr(i) = round(errMat(i,1),2);
    Ek1(i) = (activationParamsCenterDropoutEarly(i,2));
    Ek1err(i) = (errMat(i,2));
    Ek2(i) = (activationParamsCenterDropoutEarly(i,3));
    Ek2err(i) = (errMat(i,3));
end
% onTimes = table(Reporter,A,Aerr,k1,k1err,k2,k2err)
% taus:
Et1 = 1./Ek1;
Et2 = 1./Ek2;
%propagate errror:
Et1err = round(1./Ek1.^2.*Ek1err);
Et2err = round(1./Ek2.^2.*Ek2err,-1);
%slight adjustment, bc a couple of the errors are a bit crazy, and if we fit instead using taus they become reasonable:
% Et1err(1) = round(activationParamsCenterDropoutEarlyStdTaus1000(1,2));
% Et1err(3) = round(activationParamsCenterDropoutEarlyStdTaus1000(3,2),-1);
% Et2err(4) = round(activationParamsCenterDropoutEarlyStdTaus1000(4,3),-1);
% b/c we call the tau1 of neutral the tau2:
Et2(1) = Et1(1);Et2err(1) = Et1err(1);EA(1) = 0;
% round
Et1 = round(Et1);
Et2 = round(Et2,-1);
activeParamsCenterDropoutEarly = table(Reporter,EA,EAerr,Et1,Et1err,Et2,Et2err)

% 211015 updated these numbers hardenVincentDepace2021tablesV3

%% bar chart. active params. EARLY, center Fig. 3C
ax229 = figure(229);
bar(1,EA(1),0.4,'FaceColor',colr{1},'EdgeColor','none'); %initiate the bar chart first
for i = 2:5
    hold on; bar(i,EA(i),0.4,'FaceColor',colr{i},'EdgeColor','none');
end
xticks([1 2 3 4 5]);
set(gca,'xticklabels',tits(1:5),'yLim',[0 1],'Box',true,'FontSize',16,'xLim',[0.5 5.5]);
xticklabel_rotate([],45,[],'FontSize',16);
hold on;
er = errorbar(1,EA(1),EAerr(1),EAerr(1));
er.Color = colr{2};    %make the black neutral grey                       
er.LineStyle = 'none'; 
for i = 2:5
    hold on;
    er = errorbar(i,EA(i),EAerr(i),EAerr(i));
    er.Color = colr{1};    % make all err bars black                       
    er.LineStyle = 'none'; 
end
set(gca,'Box',true);
% k1's -- full distribution
ax230 = figure(230);
bar(1,Et1(1),0.4,'FaceColor',colr{1},'EdgeColor','none');
for i = 2:5
    hold on; bar(i,Et1(i),0.4,'FaceColor',colr{i},'EdgeColor','none');
end
xticks([1 2 3 4 5]);
set(gca,'xticklabels',tits(1:5),'yLim',[0 60],'Box',true,'FontSize',16,'xLim',[0.5 5.5]);
xticklabel_rotate([],45,[],'FontSize',16);
hold on;
er = errorbar(1,Et1(1),Et1err(1),Et1err(1));
er.Color = colr{2};                           
er.LineStyle = 'none'; 
for i = 2:5
    hold on;
    er = errorbar(i,Et1(i),Et1err(i),Et1err(i));
    er.Color = colr{1};                           
    er.LineStyle = 'none'; 
end
set(gca,'Box',true);
% k2's -- full distribution
ax231 = figure(231);
bar(1,Et2(1),0.4,'FaceColor',colr{1},'EdgeColor','none');
for i = 2:5
    hold on; bar(i,Et2(i),0.4,'FaceColor',colr{i},'EdgeColor','none');
end
xticks([1 2 3 4 5]);
set(gca,'xticklabels',tits(1:5),'yLim',[0 750],'Box',true,'FontSize',16,'xLim',[0.5 5.5]);
xticklabel_rotate([],45,[],'FontSize',16);
hold on;
er = errorbar(1,Et2(1),Et2err(1),Et2err(1));
er.Color = colr{2};                           
er.LineStyle = 'none'; 
for i = 2:5
    hold on;
    er = errorbar(i,Et2(i),Et2err(i),Et2err(i));
    er.Color = colr{1};                           
    er.LineStyle = 'none'; 
end
set(gca,'Box',true); 
%% save bar charts
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
% saveas(ax229,[fp 'freqBarAsEarly211015'],'svg');
% saveas(ax230,[fp 'freqBarTau1Early211015'],'svg');
% saveas(ax231,[fp 'freqBarTau2Early211015'],'svg');
saveas(ax229,[fp 'freqBarAsEarly220401'],'svg'); % include dropouts
saveas(ax230,[fp 'freqBarTau1Early220401'],'svg');
saveas(ax231,[fp 'freqBarTau2Early220401'],'svg');

%% active TRXN surv, LATE, center, (SFig. 7)
fignum = 35;
inarg = [0.9 0.01 0.001];
numOfBS = 10000;
split = 1/5;timestep = 15;
onHalf2cell1 = {};onHalf2cell2 = {};onHalf2cell3 = {};onHalf2cell4 = {};onLateParamsCell = {};
activationParamsCenterLate = zeros(5,3);
freqCenterDropoutLate = zeros(5,1);
eventNumberCenterDropoutLate = zeros(5,1);
lowTimeCenterDropoutLate = zeros(5,1);
flag = 1; 
for i = 1:5
    % treat each rep indie:
    mat1 = centerCia1{i};
    L = size(mat1,1);
    ts = sort(mat1(:,7));
    t50 = ts(floor(L*split));
    logi = mat1(:,7) > t50;
    half2mat = mat1(logi,:);
%     onHalf2cell1{i} = dropout(half2mat);
    onHalf2cell1{i} = (half2mat);
    % the other rep:
    mat2 = centerCia2{i};
    L = size(mat2,1);
    ts = sort(mat2(:,7));
    t50 = ts(floor(L*split));
    logi = mat2(:,7) > t50;
    half2mat = mat2(logi,:);
%     onHalf2cell2{i} = dropout(half2mat);
    onHalf2cell2{i} = (half2mat);
    try
        % treat each rep indie:
        mat3 = centerCia3{i};
        L = size(mat3,1);
        ts = sort(mat3(:,7));
        t50 = ts(floor(L*split));
        logi = mat3(:,7) > t50;
        half2mat = mat3(logi,:);
%         onHalf2cell3{i} = dropout(half2mat);
        onHalf2cell3{i} = (half2mat);
        % the other rep:
        mat4 = centerCia4{i};
        L = size(mat4,1);
        ts = sort(mat4(:,7));
        t50 = ts(floor(L*split));
        logi = mat4(:,7) > t50;
        half2mat = mat4(logi,:);
%         onHalf2cell4{i} = dropout(half2mat);
        onHalf2cell4{i} = (half2mat);
        % frequency:
%         [out,p] = frequency_dwellFli5(onHalf2cell1{i},onHalf2cell2{i},onHalf2cell3{i},onHalf2cell4{i},timestep,colr{i},fignum,inarg,[]);
        % for survival:
        [out,p] = survivalPlotFromVector([onHalf2cell1{i};onHalf2cell2{i};onHalf2cell3{i};onHalf2cell4{i}],6,2,1,3600,fignum,inarg,colr{i}); %
    catch
        % frequency:
%         [out,p] = frequency_dwellFli5(onHalf2cell1{i},onHalf2cell2{i},timestep,colr{i},fignum,inarg,[]);
        % for survival:
        [out,p] = survivalPlotFromVector([onHalf2cell1{i};onHalf2cell2{i}],6,2,1,3600,fignum,inarg,colr{i}); %
    end
    activationParamsCenterLate(i,:) = [out.a out.k1 out.k2];
    % frequency:
%     freqCenterDropoutLate(i) = 1/out.totalFrequency;
%     eventNumberCenterDropoutLate(i) = out.eventNumber;
%     lowTimeCenterDropoutLate(i) = out.lowTimes;
    hold on
    set(p,'Color',colr{i});
    if flag == 1
        try
            % frequency:
%             [~,pbs] = bootstrapsFrequency([onHalf2cell1{i};onHalf2cell2{i};onHalf2cell3{i};onHalf2cell4{i}],6,numOfBS,freqCenterDropoutLate(i),colr{i},fignum);
            % survival:
            [~,pbs] = bootstrapsBootstraps([onHalf2cell1{i};onHalf2cell2{i};onHalf2cell3{i};onHalf2cell4{i}],6,numOfBS,colr{i},fignum);
        catch
            % frequency:
%             [~,pbs] = bootstrapsFrequency([onHalf2cell1{i};onHalf2cell2{i}],6,numOfBS,freqCenterDropoutLate(i),colr{i},fignum);
            % survial:
            [~,pbs] = bootstrapsBootstraps([onHalf2cell1{i};onHalf2cell2{i}],6,numOfBS,colr{i},fignum);
        end
        hold on
    end
end
% frequency:
% set(gca,'YLim',[1.0e-5 0.02],'XLim',[0 1700]); % all cndns
% for survival:
set(gca,'Box',true,'FontSize',16,'YLim',[-6.05 0],'XLim',[0 1700]);
%% save
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
% saveas(gca,[fp 'survFreqCurveFit&errLate220328'],'svg'); %saved with 10000 BSs
% saveas(gca,[fp 'survFreqCurveFit&errLate220401'],'svg'); % frequency with dropouts
saveas(gca,[fp 'survFreqCurveFit&errLate220413'],'svg'); % survival curves
%% error surv Freq, LATE, center, SFIg. 7
fitParams = [1 3600 0.9 0.01 0.001]; %[tm tx a k1 k2]
bsNum = 10000;
% parpool
tic;
activationParamsCenterDropoutLateStd1000 = [];activeCenterLateFits = {};activeCenterLateFitsTaus = {};
activationParamsCenterDropoutLateStdTaus1000 = [];mat ={};
% parfor i = 1:5
%     size(non_crack_bytes); 
for i = 1:5
    try
        mat{i} = [onHalf2cell1{i}; onHalf2cell2{i}; onHalf2cell3{i}; onHalf2cell4{i}];
    catch
        mat{i} = [onHalf2cell1{i}; onHalf2cell2{i}];
    end
    activeCenterLateFitsTaus{i} = bootstrapsBootstrapsFitTaus(mat{i},6,bsNum,2,fitParams);
    activationParamsCenterDropoutLateStdTaus1000(i,:) = activeCenterLateFitsTaus{i}.std;
    activeCenterLateFits{i} = bootstrapsBootstrapsFit(mat{i},6,bsNum,2,fitParams);
    activationParamsCenterDropoutLateStd1000(i,:) = activeCenterLateFits{i}.std;
end
toc

%% save the Bootstraps:
save C:\Users\thard\Dropbox\matlab\data\freqSurvLateErrs220328.dat activeCenterLateFitsTaus activationParamsCenterDropoutLateStdTaus1000 activeCenterLateFits activationParamsCenterDropoutLateStd1000

%% table. active params. LATE, center, TABLE 2
% set up table
Reporter = tits';
% rates:
% allocate:
A = zeros(5,1);k1 = A; k2 = A; %allocate space, amke things col vectors
Aerr = A; k1err = A; k2err = A;
errMat = activationParamsCenterDropoutLateStd1000;
for i = 1:5 
    A(i) = round(activationParamsCenterDropoutLate(i,1),2);
    Aerr(i) = round(errMat(i,1),2);
    k1(i) = (activationParamsCenterDropoutLate(i,2));
    k1err(i) = (errMat(i,2));
    k2(i) = (activationParamsCenterDropoutLate(i,3));
    k2err(i) = (errMat(i,3));
end
% onTimes = table(Reporter,A,Aerr,k1,k1err,k2,k2err)
% taus:
t1 = 1./k1;
t2 = 1./k2;
%propagate errror:
t1err = round(1./k1.^2.*k1err);
t2err = round(1./k2.^2.*k2err,-1);
% round
t1 = round(t1);
t2 = round(t2,-1);
activeParamsCenterDropoutLate = table(Reporter,A,Aerr,t1,t1err,t2,t2err)
% saved to hardenVincentDepace2021tablesV3 on 211018

%% comparing replicates -- spots v time (unused)
%Importantly, here we normalize by the median number of nukes in all
%experiments. Run the nukes binning section for that
% neutral:
ax30 = figure(30);
for j = 1:4 % bc 4 reps
    %get all rose with a spot:
    mat0 = nNukesCell{j};
    mat = mat0(mat0(:,4) > 0,:); %col 4 is binary
    size(unique(mat(:,3)));
    %get frame range:
    mx = max(mat(:,2)); %frame indexing is more useful that time
    svt = []; %spots v time
    for i = 1:mx
        mat2 = mat(mat(:,2) == i,:); %find all rose for this frame (logical indexing)
        if ~isempty(mat2) %if we have nukes with spots, mark the time and sum up the num of spots
            svt = [svt;mat2(1,1) size(mat2,1)]; %[time numOfSpots]
        end
    end
    hold on;plot(svt(:,1),svt(:,2)./(medNukeNum/2),'Color',colr{1} + j*0.2,'LineStyle','-','Marker','.','MarkerSize',10);shg
end
set(gca,'Box',true,'FontSize',16,'xLim',[0 4000]);

%% wt space:
ax31 = figure(31);
spaceSVT = {};
for j = 1:2 % bc 4 reps
    %get all rose with a spot:
    mat0 = spaceNukesCell{j};
    mat = mat0(mat0(:,4) > 0,:); %col 4 is binary
    size(unique(mat(:,3)));
    %get frame range:
    mx = max(mat(:,2)); %frame indexing is more useful that time
    svt = []; %spots v time
    for i = 1:mx
        mat2 = mat(mat(:,2) == i,:); %find all rose for this frame (logical indexing)
        if ~isempty(mat2) %if we have nukes with spots, mark the time and sum up the num of spots
            svt = [svt;mat2(1,1) size(mat2,1)]; %[time numOfSpots]
        end
    end
    hold on;plot(svt(:,1),svt(:,2)./(medNukeNum/2),'Color',colr{1} + (j - 1)*0.5,'LineStyle','-','Marker','.','MarkerSize',10);shg
    spaceSVT{j} = svt(:,2);
end
set(gca,'Box',true,'FontSize',16,'xLim',[0 4000]);
[h,p] = kstest2(spaceSVT{1},spaceSVT{2}) % these attemps at a KS test are largly a failure, as they fail for all but Bcd. I am not even sure this is the correct way to use ks tests, as these are time dependent vectors, not collections of intervals or whatever

%% zld:
ax32 = figure(32);
zldSVT = {};
for j = 1:2 % bc 4 reps
    %get all rose with a spot:
    mat0 = zldNukesCell{j};
    mat = mat0(mat0(:,4) > 0,:); %col 4 is binary
    size(unique(mat(:,3)));
    %get frame range:
    mx = max(mat(:,2)); %frame indexing is more useful that time
    svt = []; %spots v time
    for i = 1:mx
        mat2 = mat(mat(:,2) == i,:); %find all rose for this frame (logical indexing)
        if ~isempty(mat2) %if we have nukes with spots, mark the time and sum up the num of spots
            svt = [svt;mat2(1,1) size(mat2,1)]; %[time numOfSpots]
        end
    end
    hold on;plot(svt(:,1),svt(:,2)./(medNukeNum/2),'Color',colr{1} + (j - 1)*0.5,'LineStyle','-','Marker','.','MarkerSize',10);shg
    zldSVT{j} = svt(:,2);
end
set(gca,'Box',true,'FontSize',16,'xLim',[0 4000]);
[h,p] = kstest2(zldSVT{1},zldSVT{2})
%% bcd:
ax33 = figure(33);
bcdSVT = {};
for j = 1:2 % bc 4 reps
    %get all rose with a spot:
    mat0 = bcdNukesCell{j};
    mat = mat0(mat0(:,4) > 0,:); %col 4 is binary
    size(unique(mat(:,3)));
    %get frame range:
    mx = max(mat(:,2)); %frame indexing is more useful that time
    svt = []; %spots v time
    for i = 1:mx
        mat2 = mat(mat(:,2) == i,:); %find all rose for this frame (logical indexing)
        if ~isempty(mat2) %if we have nukes with spots, mark the time and sum up the num of spots
            svt = [svt;mat2(1,1) size(mat2,1)]; %[time numOfSpots]
        end
    end
    hold on;plot(svt(:,1),svt(:,2)./(medNukeNum/2),'Color',colr{1} + (j - 1)*0.5,'LineStyle','-','Marker','.','MarkerSize',10);shg
    bcdSVT{j} = svt(:,2);
end
set(gca,'Box',true,'FontSize',16,'xLim',[0 4000]);
[h,p] = kstest2(bcdSVT{1},bcdSVT{2})
%% d:
ax34 = figure(34);
for j = 1:2 % bc 4 reps
    %get all rose with a spot:
    mat0 = dNukesCell{j};
    mat = mat0(mat0(:,4) > 0,:); %col 4 is binary
    size(unique(mat(:,3)));
    %get frame range:
    mx = max(mat(:,2)); %frame indexing is more useful that time
    svt = []; %spots v time
    for i = 1:mx
        mat2 = mat(mat(:,2) == i,:); %find all rose for this frame (logical indexing)
        if ~isempty(mat2) %if we have nukes with spots, mark the time and sum up the num of spots
            svt = [svt;mat2(1,1) size(mat2,1)]; %[time numOfSpots]
        end
    end
    hold on;plot(svt(:,1),svt(:,2)./(medNukeNum/2),'Color',colr{1} + (j - 1)*0.5,'LineStyle','-','Marker','.','MarkerSize',10);shg
end
set(gca,'Box',true,'FontSize',16,'xLim',[0 4000]);

%% save
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(ax30,[fp 'trxnProfileReplicatesNeutral211108'],'svg');
saveas(ax31,[fp 'trxnProfileReplicatesWT211108'],'svg');
saveas(ax32,[fp 'trxnProfileReplicatesZld211108'],'svg');
saveas(ax33,[fp 'trxnProfileReplicatesBcd211108'],'svg');
saveas(ax34,[fp 'trxnProfileReplicatesDst211108'],'svg');
% saved 211108

% On biological replicates:

%% KS test btwn replicates for on frames, center nukes:
% order: neutral, wt space, zld, bcd, stat
for i = 2:5
    m1 = dropout(centerCia1{i});m2 = dropout(centerCia2{i});
    tits{i} % to make the output easier to read
    [h,p] = kstest2(m1(:,5),m2(:,5))
end
% fails to reject the null hypothesis for all conditions at p = 0.05

%% KS test for first passage frames. Central nukes:
for i = 2:5
    tits{i} % to make the output easier to read
    [v1,v2] = firstPassageReplicates(analCell1{i},analCell2{i},bins,boiC);
%     [h,p] = kstest2(v1,v2,'alpha',0.01)
    [h,p] = kstest2(v1,v2)
%     avg = mean([size(v1,1) size(v2,1)])  % report the number of nuclei per construct if you want
%     stDev = std([size(v1,1) size(v2,1)])
end
% order: neutral, wt space, zld, bcd, stat

% fails to rejecrt the null hypothesis for all conditions but neutral

%% KS test for Idle times
for i = 2:5
    interEvents1 = idleSurvival(centerCia1{i});
    interEvents2 = idleSurvival(centerCia2{i});
    [h,p] = kstest2(interEvents1,interEvents2)
end
% fails to reject the null hypothesis for all conditions at p = 0.05

% KS test for many neutral replicates:
%% first passages KS for all 4 neutral reps:
% first get the vectors of events:
[v1,v2] = firstPassageReplicates(analCell1{1},analCell2{1},bins,boiC);
[v3,v4] = firstPassageReplicates(analCell3{1},analCell4{1},bins,boiC);
%     [h,p] = kstest2(v1,v2,'alpha',0.01)
% now compare each replicate against the collection of the others:
[h,p] = kstest2(v1,[v2;v3;v4])
[h,p] = kstest2(v2,[v1;v3;v4]) % fails
[h,p] = kstest2(v3,[v1;v2;v4])
[h,p] = kstest2(v4,[v1;v2;v3]) % fails

%% how many active nuclei in each neutral embryo?
mean([size(v1,1) size(v2,1) size(v3,1) size(v4,1)])

%% active trxn KS for all 4 neutral reps:
% get the distributions
m1 = dropout(centerCia1{1});m2 = dropout(centerCia2{1});m3 = dropout(centerCia3{1});m4 = dropout(centerCia4{1});
% now compare each replicate against the collection of the others:
[h,p] = kstest2(m1(:,5),[m2(:,5);m3(:,5);m4(:,5)])
[h,p] = kstest2(m2(:,5),[m1(:,5);m3(:,5);m4(:,5)],'alpha',0.01) % fails at p = 0.05
[h,p] = kstest2(m3(:,5),[m1(:,5);m2(:,5);m4(:,5)])
[h,p] = kstest2(m4(:,5),[m1(:,5);m2(:,5);m3(:,5)]) 

%% first passage time by replicate (SFig. 4A):
% neutral:
% plot data first:
fignum = 117;
ax117 = figure(fignum);
repNum = 4;
analCell = nAnalCell;
ciaCell = nCenterCia;
NtV = zeros(1,repNum);
for i = 1:repNum
    [out,p] = initialFractionBin2(analCell{i},bins,boiC,[],30,3000,fignum);
    set(p,'Color',colr{1} + i*0.2,'LineWidth',1);
    hold on
    NtV(i) = out.totNukes;
end
set(gca,'Box',true,'FontSize',16,'xLim',[0 4000],'yLim',[0 1.01]);
% fit:
inarg = [200 12 12 12 12 0.4 0.4 0.4 0.4];
% inarg = [200 12 0.4]; % inarg = [200 12 4 4 8 8 0.4    0.9    0.9    0.7    0.9]; 
nFPrepFit = fminsearch('globalGammaFit5',inarg,[],nCenterCia,NtV);
% add fits to plot:
xs =0:4000;
rate = nFPrepFit(1);
stepsV = nFPrepFit(2:repNum + 1);
AfV = nFPrepFit(repNum + 2:end);
for i = 1:repNum
    gCDF = AfV(i)*gamcdf(xs,stepsV(i),rate);
    hold on;plot(xs,gCDF,'--','LineWidth',0.5,'Color',colr{1} + i*0.2,'LineWidth',1);shg
end
%% save
fp = 'C:\Users\thard\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(ax117,[fp 'fpReplicatesNeutral211111'],'svg');
%% neutral rep first bind BS
% generate bootstraps
bs = 10000; 
nBsCell = {};
tic;
for i = 1:4 % for 4 replicates
    nBsCell{i} = firstEventBootstrapGenerator2(nAnalCell{i},bins,boiC,bs);
end
toc
% nBS fit:
inarg = [200 12 12 12 12 0.4 0.4 0.4 0.4];
tic;
nFitErrStruc = firstEventBootstrapFitter3(nBsCell,inarg);
toc

%% saved so you dont ahve to run again:
% save C:\Users\thard\Dropbox\matlab\data\t1p54n.dat nFitErrStruc nBsCell
save C:\Users\thard\Dropbox\matlab\data\t1p54n.dat nFitErrStruc nBsCell -append
%% plot neutral bootstraps
fignum = 118; 
ax118 = figure(fignum);
% load C:\Users\thard\Dropbox\matlab\data\t1p54n1.dat -mat 
tic;
for i = 1:4
    [~,p] = firstEventBootstrapPlotter2(nBsCell{i},nAnalCell{i},colr{1} + i*0.2,fignum);
    hold on;
end
set(gca,'YLim',[0 1.01],'Xlim',[0 4000],'Box',true,'FontSize',16); 
% ylabel('Fraction of nuclei')
% xlabel('Time (s)')
toc;
%% save 
fp = 'C:\Users\thard\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(ax118,[fp 'fpReplicatesNeutralbs211111'],'svg');
% done 211112

%% space:
fignum = 119;
ax119 = figure(fignum);
analCell = spaceAnalCell;
ciaCell = spaceCenterCia;
repNum = 2;
NtV = zeros(1,repNum);
for i = 1:repNum
    [out,p] = initialFractionBin2(analCell{i},bins,boiC,[],30,3000,fignum);
    set(p,'Color',colr{1} + (i - 1)*0.5,'LineWidth',1);
    hold on
    NtV(i) = out.totNukes;
end
set(gca,'Box',true,'FontSize',16,'xLim',[0 4000],'yLim',[0 1]);
% fit:
inarg = [200 4 4 0.9 0.9];
spaceFPrepFit = fminsearch('globalGammaFit5',inarg,[],ciaCell,NtV);
% add fits to plot:
xs =0:4000;
rate = spaceFPrepFit(1);
stepsV = spaceFPrepFit(2:repNum + 1);
AfV = spaceFPrepFit(repNum + 2:end);
for i = 1:repNum
    gCDF = AfV(i)*gamcdf(xs,stepsV(i),rate);
    hold on;plot(xs,gCDF,'--','LineWidth',0.5,'Color',colr{1} + (i - 1)*0.5,'LineWidth',1);shg
end
%% save
fp = 'C:\Users\thard\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(ax119,[fp 'fpReplicatesWTspace211111'],'svg');
%% space rep first bind BS
% generate bootstraps
bs = 10000; 
spaceBsCell = {};
tic;
for i = 1:2 % for 2 replicates
    spaceBsCell{i} = firstEventBootstrapGenerator2(spaceAnalCell{i},bins,boiC,bs);
end
toc
% space BS fit:
inarg = [200 4 4 0.9 0.9];
tic;
spaceFitErrStruc = firstEventBootstrapFitter3(spaceBsCell,inarg);
toc

%% saved so you dont ahve to run again:
save C:\Users\thard\Dropbox\matlab\data\t1p54n.dat spaceFitErrStruc spaceBsCell -append
%%
% plot space bootstraps
fignum = 120; 
ax120 = figure(fignum);
% load C:\Users\thard\Dropbox\matlab\data\t1p54n.dat -mat 
tic;
for i = 1:2
    [~,p] = firstEventBootstrapPlotter2(spaceBsCell{i},spaceAnalCell{i},colr{1} + (i - 1)*0.5,fignum);
    hold on;
end
set(gca,'YLim',[0 1.01],'Xlim',[0 4000],'Box',true,'FontSize',16); 
% ylabel('Fraction of nuclei')
% xlabel('Time (s)')
toc;
%% save 
fp = 'C:\Users\thard\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(ax120,[fp 'fpReplicatesWTspacebs211111'],'svg');
% done 211112

%% zld:
fignum = 121;
ax121 = figure(fignum);
analCell = zldAnalCell;
ciaCell = zldCenterCia;
repNum = 2;
NtV = zeros(1,repNum);
for i = 1:repNum
    [out,p] = initialFractionBin2(analCell{i},bins,boiC,[],30,3000,fignum);
    set(p,'Color',colr{1} + (i - 1)*0.5,'LineWidth',1);
    hold on
    NtV(i) = out.totNukes;
end
set(gca,'Box',true,'FontSize',16,'xLim',[0 4000],'yLim',[0 1]);
% fit:
inarg = [200 4 4 0.9 0.9];
zldFPrepFit = fminsearch('globalGammaFit5',inarg,[],ciaCell,NtV);
% add fits to plot:
xs =0:4000;
rate = zldFPrepFit(1);
stepsV = zldFPrepFit(2:repNum + 1);
AfV = zldFPrepFit(repNum + 2:end);
for i = 1:repNum
    gCDF = AfV(i)*gamcdf(xs,stepsV(i),rate);
    hold on;plot(xs,gCDF,'--','LineWidth',0.5,'Color',colr{1} + (i - 1)*0.5,'LineWidth',1);shg
end
%% save
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(ax121,[fp 'fpReplicatesZld220329'],'svg');
%% zld rep first bind BS
% generate bootstraps
bs = 10000; 
zldBsCell = {};
tic;
for i = 1:2 % for 2 replicates
    zldBsCell{i} = firstEventBootstrapGenerator2(zldAnalCell{i},bins,boiC,bs);
end
toc
% zld BS fit:
inarg = [200 4 4 0.9 0.9];
tic;
zldFitErrStruc = firstEventBootstrapFitter3(zldBsCell,inarg);
toc

%% saved so you dont ahve to run again:
save C:\Users\tth12\Dropbox\matlab\data\t1p54n.dat zldFitErrStruc zldBsCell -append
%%
% plot zld bootstraps
fignum = 122; 
ax122 = figure(fignum);
% load C:\Users\thard\Dropbox\matlab\data\t1p54n.dat -mat 
tic;
for i = 1:2
    [~,p] = firstEventBootstrapPlotter2(zldBsCell{i},zldAnalCell{i},colr{1} + (i - 1)*0.5,fignum);
    hold on;
end
set(gca,'YLim',[0 1.01],'Xlim',[0 4000],'Box',true,'FontSize',16); 
% ylabel('Fraction of nuclei')
% xlabel('Time (s)')
toc;
%% save 
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(ax122,[fp 'fpReplicatesZldBS220329'],'svg');
% done 211112

%% bcd:
fignum = 123;
ax123 = figure(fignum);
analCell = bcdAnalCell;
ciaCell = bcdCenterCia;
repNum = 2;
NtV = zeros(1,repNum);
for i = 1:repNum
    [out,p] = initialFractionBin2(analCell{i},bins,boiC,[],30,3000,fignum);
    set(p,'Color',colr{1} + (i - 1)*0.5,'LineWidth',1);
    hold on
    NtV(i) = out.totNukes;
end
set(gca,'Box',true,'FontSize',16,'xLim',[0 4000],'yLim',[0 1]);
% fit:
inarg = [200 8 8 0.7 0.7];
bcdFPrepFit = fminsearch('globalGammaFit5',inarg,[],ciaCell,NtV);
% add fits to plot:
xs =0:4000;
rate = bcdFPrepFit(1);
stepsV = bcdFPrepFit(2:repNum + 1);
AfV = bcdFPrepFit(repNum + 2:end);
for i = 1:repNum
    gCDF = AfV(i)*gamcdf(xs,stepsV(i),rate);
    hold on;plot(xs,gCDF,'--','LineWidth',0.5,'Color',colr{1} + (i - 1)*0.5,'LineWidth',1);shg
end
%% save
fp = 'C:\Users\thard\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(ax123,[fp 'fpReplicatesBcd211111'],'svg');
%% bcd rep first bind BS
% generate bootstraps
bs = 10000; 
bcdBsCell = {};
tic;
for i = 1:2 % for 2 replicates
    bcdBsCell{i} = firstEventBootstrapGenerator2(bcdAnalCell{i},bins,boiC,bs);
end
toc
% bcd BS fit:
inarg = [200 8 8 0.7 0.7];
tic;
bcdFitErrStruc = firstEventBootstrapFitter3(bcdBsCell,inarg);
toc

%% saved so you dont ahve to run again:
save C:\Users\thard\Dropbox\matlab\data\t1p54n.dat bcdFitErrStruc bcdBsCell -append
%%
% plot bcd bootstraps
fignum = 124; 
ax124 = figure(fignum);
% load C:\Users\thard\Dropbox\matlab\data\t1p54n.dat -mat 
tic;
for i = 1:2
    [~,p] = firstEventBootstrapPlotter2(bcdBsCell{i},bcdAnalCell{i},colr{1} + (i - 1)*0.5,fignum);
    hold on;
end
set(gca,'YLim',[0 1.01],'Xlim',[0 4000],'Box',true,'FontSize',16); 
% ylabel('Fraction of nuclei')
% xlabel('Time (s)')
toc;
%% save 
fp = 'C:\Users\thard\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(ax124,[fp 'fpReplicatesbcdBS211111'],'svg');

%% stat:
fignum = 125;
ax125 = figure(fignum);
analCell = dAnalCell;
ciaCell = dCenterCia;
repNum = 2;
NtV = zeros(1,repNum);
for i = 1:repNum
    [out,p] = initialFractionBin2(analCell{i},bins,boiC,[],30,3000,fignum);
    set(p,'Color',colr{1} + (i - 1)*0.5,'LineWidth',1);
    hold on
    NtV(i) = out.totNukes;
end
set(gca,'Box',true,'FontSize',16,'xLim',[0 4000],'yLim',[0 1]);
%TS ing 
NtV = [50 45]; % not sure why, but the fitting algorithm didnt like when we got really close to Af = 1. So I bumped up the number of nukes by 1 for each rep and we got well behaved fitting. 
% fit:
inarg = [200 8 8 0.9 0.9];
dFPrepFit = fminsearch('globalGammaFit5',inarg,[],ciaCell,NtV);
% add fits to plot:
xs =0:4000;
rate = dFPrepFit(1);
stepsV = dFPrepFit(2:repNum + 1);
AfV = dFPrepFit(repNum + 2:end);
for i = 1:repNum
    gCDF = AfV(i)*gamcdf(xs,stepsV(i),rate);
    hold on;plot(xs,gCDF,'--','LineWidth',0.5,'Color',colr{1} + (i - 1)*0.5,'LineWidth',1);shg
end
%% save
fp = 'C:\Users\thard\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(ax125,[fp 'fpReplicatesStat211111'],'svg');
%% stat rep first bind BS
% generate bootstraps
bs = 10000; 
dBsCell = {};
tic;
for i = 1:2 % for 2 replicates
    dBsCell{i} = firstEventBootstrapGenerator2(dAnalCell{i},bins,boiC,bs);
end
toc
% stat BS fit:
inarg = [200 8 8 0.9 0.9];
tic;
dFitErrStruc = firstEventBootstrapFitter3(dBsCell,inarg);
toc

%% saved so you dont ahve to run again:
save C:\Users\thard\Dropbox\matlab\data\t1p54n.dat dFitErrStruc dBsCell -append
%%
% plot d bootstraps
fignum = 126; 
ax126 = figure(fignum);
% load C:\Users\thard\Dropbox\matlab\data\t1p54n.dat -mat 
tic;
for i = 1:2
    [~,p] = firstEventBootstrapPlotter2(dBsCell{i},dAnalCell{i},colr{1} + (i - 1)*0.5,fignum);
    hold on;
end
set(gca,'YLim',[0 1.01],'Xlim',[0 4000],'Box',true,'FontSize',16); 
% ylabel('Fraction of nuclei')
% xlabel('Time (s)')
toc;
%% save 
fp = 'C:\Users\thard\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(ax126,[fp 'fpReplicatesdBS211111'],'svg');
% done 211112 & ready to make this part of the Sfigure. 

%% active trxn time by replicate (SFig. 4):
%  neutral: with the fit parameters
fignum = 130;
ciaCell = nCenterCia;
repNum = 4;
numBS = 10000;
flag = 1;
for i = 1:repNum
%     [~,p] = survivalPlotFromVector(dropout(ciaCell{i}),6,0,1,3600,fignum + 1); % no fit
    [~,p] = survivalPlotFromVector((ciaCell{i}),6,0,1,3600,fignum + 1); % no fit
    set(p,'Color',colr{1} + i*0.2);
    hold on
    if flag == 1
%         bootstrapsBootstraps(dropout(ciaCell{i}),6,numBS,colr{1} + i*0.2,fignum + 1);
        bootstrapsBootstraps((ciaCell{i}),6,numBS,colr{1} + i*0.2,fignum + 1);
        hold on
    end
end
set(gca,'Box',true,'FontSize',16,'YLim',[-6 0],'XLim',[0 2500]);
% -4.3, 600
% space:  with the fit parameters
fignum = 131;
ciaCell = spaceCenterCia;
repNum = 2;
numBS = 10000;
flag = 1;
for i = 1:repNum
%     [~,p] = survivalPlotFromVector(dropout(ciaCell{i}),6,0,1,3600,fignum + 1); % no fit
    [~,p] = survivalPlotFromVector((ciaCell{i}),6,0,1,3600,fignum + 1); % no fit
    set(p,'Color',colr{i});
    hold on
    if flag == 1
%         bootstrapsBootstraps(dropout(ciaCell{i}),6,numBS,colr{i},fignum + 1);
        bootstrapsBootstraps((ciaCell{i}),6,numBS,colr{i},fignum + 1);
        hold on
    end
end
set(gca,'Box',true,'FontSize',16,'YLim',[-6 0],'XLim',[0 2500]);
% -6, 2500
% zld:
fignum = 132;
ciaCell = zldCenterCia;
repNum = 2;
numBS = 10000;
flag = 1;
for i = 1:repNum
%     [~,p] = survivalPlotFromVector(dropout(ciaCell{i}),6,0,1,3600,fignum + 1); % no fit
    [~,p] = survivalPlotFromVector((ciaCell{i}),6,0,1,3600,fignum + 1); % no fit
    set(p,'Color',colr{i});
    hold on
    if flag == 1
%         bootstrapsBootstraps(dropout(ciaCell{i}),6,numBS,colr{i},fignum + 1);
        bootstrapsBootstraps((ciaCell{i}),6,numBS,colr{i},fignum + 1);
        hold on
    end
end
set(gca,'Box',true,'FontSize',16,'YLim',[-6 0],'XLim',[0 2500]);
% -5.3; 1800
% bcd:
fignum = 133;
ciaCell = bcdCenterCia;
repNum = 2;
numBS = 10000;
flag = 1;
for i = 1:repNum
%     [~,p] = survivalPlotFromVector(dropout(ciaCell{i}),6,0,1,3600,fignum + 1); % no fit
    [~,p] = survivalPlotFromVector((ciaCell{i}),6,0,1,3600,fignum + 1); % no fit
    set(p,'Color',colr{i});
    hold on
    if flag == 1
%         bootstrapsBootstraps(dropout(ciaCell{i}),6,numBS,colr{i},fignum + 1);
        bootstrapsBootstraps((ciaCell{i}),6,numBS,colr{i},fignum + 1);
        hold on
    end
end
set(gca,'Box',true,'FontSize',16,'YLim',[-6 0],'XLim',[0 2500]);
% -5; 900
% dstat:
fignum = 134;
ciaCell = dCenterCia;
repNum = 2;
numBS = 10000;
flag = 1;
for i = 1:repNum
%     [~,p] = survivalPlotFromVector(dropout(ciaCell{i}),6,0,1,3600,fignum + 1); % no fit
    [~,p] = survivalPlotFromVector((ciaCell{i}),6,0,1,3600,fignum + 1); % no fit
    set(p,'Color',colr{i});
    hold on
    if flag == 1
%         bootstrapsBootstraps(dropout(ciaCell{i}),6,numBS,colr{i},fignum + 1);
        bootstrapsBootstraps((ciaCell{i}),6,numBS,colr{i},fignum + 1);
        hold on
    end
end
set(gca,'Box',true,'FontSize',16,'YLim',[-6 0],'XLim',[0 2500]);
% -5.5; 2200
%% save
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
figure(131);saveas(gca,[fp 'survCurveFit&errNeutralReps220502'],'svg');
figure(132);saveas(gca,[fp 'survCurveFit&errWtReps220502'],'svg');
figure(133);saveas(gca,[fp 'survCurveFit&errZldReps220502'],'svg');
figure(134);saveas(gca,[fp 'survCurveFit&errBcdReps220502'],'svg');
figure(135);saveas(gca,[fp 'survCurveFit&errDstReps220502'],'svg');
%%  neutral with bootstraps 
timestep = 15;
fignum = 136;
reps = 1:4;
ciaCell = nCenterCia;
bsNum = 10000;
% generate and plot the boot straps
tic;
for i = reps
    [~,pbs] = bootstrapsFrequency(dropout(ciaCell{i}),6,bsNum,freqCenterDropoutRep(cndnN,i),colr{1} + i*0.2,fignum);
    hold on
end
toc
% add the curves and fits:
for i = reps
    out = frequency_dwellFli3(dropout(ciaCell{i}),timestep,colr{1} + i*0.2,fignum,[0.9 0.01 0.001],[]);
    hold on
end
set(gca,'YLim',[7e-6 5e-3],'XLim',[0 2300]);
%% save this boot strap plot:
fp = 'C:\Users\thard\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(gca,[fp 'survFreqCurveFit&errNeutralReps211112'],'svg');
% done 211112

%%  space with bootstraps
timestep = 15;
fignum = 137;
reps = 1:2;
ciaCell = spaceCenterCia;
cndnN = 2;
bsNum = 10000;
% generate and plot the boot straps
tic;
for i = reps
    [~,pbs] = bootstrapsFrequency(dropout(ciaCell{i}),6,bsNum,freqCenterDropoutRep(cndnN,i),colr{1} + (i - 1)*0.5,fignum);
    hold on
end
toc
% add the curves and fits:
for i = reps
    out = frequency_dwellFli3(dropout(ciaCell{i}),timestep,colr{1} + (i - 1)*0.5,fignum,[0.9 0.01 0.001],[]);
    hold on
end
set(gca,'YLim',[7e-6 5e-3],'XLim',[0 2300]);
%% save this boot strap plot:
fp = 'C:\Users\thard\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(gca,[fp 'survFreqCurveFit&errWTspaceReps211112'],'svg');
% done 211112

%%  zld with bootstraps
timestep = 15;
fignum = 138;
reps = 1:2;
ciaCell = zldCenterCia;
cndnN = 3;
bsNum = 10000;
% generate and plot the boot straps
tic;
for i = reps
    [~,pbs] = bootstrapsFrequency(dropout(ciaCell{i}),6,bsNum,freqCenterDropoutRep(cndnN,i),colr{1} + (i - 1)*0.5,fignum);
    hold on
end
toc
% add the curves and fits:
for i = reps
    out = frequency_dwellFli3(dropout(ciaCell{i}),timestep,colr{1} + (i - 1)*0.5,fignum,[0.9 0.01 0.001],[]);
    hold on
end
set(gca,'YLim',[7e-6 5e-3],'XLim',[0 2300]);
%% save this boot strap plot:
fp = 'C:\Users\thard\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(gca,[fp 'survFreqCurveFit&errZldReps211112'],'svg');
% done 211112

%%  bcd with bootstraps
timestep = 15;
fignum = 139;
reps = 1:2;
ciaCell = bcdCenterCia;
cndnN = 4;
bsNum = 10000;
% generate and plot the boot straps
tic;
for i = reps
    [~,pbs] = bootstrapsFrequency(dropout(ciaCell{i}),6,bsNum,freqCenterDropoutRep(cndnN,i),colr{1} + (i - 1)*0.5,fignum);
    hold on
end
toc
% add the curves and fits:
for i = reps
    out = frequency_dwellFli3(dropout(ciaCell{i}),timestep,colr{1} + (i - 1)*0.5,fignum,[0.9 0.01 0.001],[]);
    hold on
end
set(gca,'YLim',[7e-6 5e-3],'XLim',[0 2300]);
%% save this boot strap plot:
fp = 'C:\Users\thard\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(gca,[fp 'survFreqCurveFit&errBcdReps211112'],'svg');
% done 211112

%%  dStat with bootstraps
timestep = 15;
fignum = 140;
reps = 1:2;
ciaCell = dCenterCia;
cndnN = 5;
bsNum = 10000;
% generate and plot the boot straps
tic;
for i = reps
    [~,pbs] = bootstrapsFrequency(dropout(ciaCell{i}),6,bsNum,freqCenterDropoutRep(cndnN,i),colr{1} + (i - 1)*0.5,fignum);
    hold on
end
toc
% add the curves and fits:
for i = reps
    out = frequency_dwellFli3(dropout(ciaCell{i}),timestep,colr{1} + (i - 1)*0.5,fignum,[0.9 0.01 0.001],[]);
    hold on
end
set(gca,'YLim',[7e-6 5e-3],'XLim',[0 2300]);
%% save this boot strap plot:
fp = 'C:\Users\thard\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(gca,[fp 'survFreqCurveFit&errStatReps211112'],'svg');
% done 211112

%NB (211112): We have all the FP & active TRXN time plots with fits &
%bootstraps. Now we just need to assemble them into a figure, replace the
%current SFig on replicates, and write something in the results and
%discussion as per the responses to reviewers

% I have no idea why this works to create a manipuble vector figure and the
% method for doing the same in FP times does not. 
%% plotting the stripe center. SFig. 3
% you could normalize by totalBinMean from the nuke binning section at the
% top of this script
lim = [-0.2 0.25];
xName = 'Fraction of nuclei';

% neutral 
figNum = 20;
anal1 = nAnal1;
anal2 = nAnal2;
anal3 = nAnal3;
anal4 = nAnal4;
col = colr{1};
ax20 = figure(figNum);
% rep 1
out1 = positionDistribution(anal1,bins);
dist1 = out1.eventBin./(out1.totalBin);dist1(isnan(dist1)) = 0; 
% rep 2
out2 = positionDistribution(anal2,bins);
dist2 = out2.eventBin./(out2.totalBin);dist2(isnan(dist2)) = 0; 
% rep 3 
out3 = positionDistribution(anal3,bins);
dist3 = out3.eventBin./(out3.totalBin);dist3(isnan(dist3)) = 0; 
% rep 4
out4 = positionDistribution(anal4,bins);
dist4 = out4.eventBin./(out4.totalBin);dist4(isnan(dist4)) = 0; 
% combined 
% handle = bar(bins,mean([dist1 dist2],2),'histc');shg
% handle = bar(bins,mean([dist1 dist2 dist3],2),'histc');shg
handle = bar(bins,mean([dist1 dist2 dist3 dist4],2),'histc');shg
patch('Parent',gca,'Vertices',handle.Vertices,'Faces',handle.Faces,...
    'FaceColor',col,...
    'EdgeColor',colr{2},...
    'CData',handle.CData);
set(gca,'Box',true,'FontSize',16,'XLim',lim,'YLim',[0 1]); 
ylabel(xName)
xlabel('Position relative to stripe center (fraction of AP length)')

%% space
figNum = 21;
anal1 = spaceAnal1;
anal2 = spaceAnal2;
col = colr{2};
ax21 = figure(figNum);
% rep 1
out1 = positionDistribution(anal1,bins);
dist1 = out1.eventBin./(out1.totalBin);dist1(isnan(dist1)) = 0; 
% rep 2
out2 = positionDistribution(anal2,bins);
dist2 = out2.eventBin./(out2.totalBin);dist2(isnan(dist2)) = 0; 
% combined 
handle = bar(bins,mean([dist1 dist2],2),'histc');shg
patch('Parent',gca,'Vertices',handle.Vertices,'Faces',handle.Faces,...
    'FaceColor',col,...
    'CData',handle.CData);
set(gca,'Box',true,'FontSize',16,'XLim',lim); 
ylabel(xName)
xlabel('Position relative to stripe center (fraction of AP length)')

%% zld
figNum = 22;
anal1 = zldAnal1;
anal2 = zldAnal2;
col = colr{3};
ax22 = figure(figNum);
% rep 1
out1 = positionDistribution(anal1,bins);
dist1 = out1.eventBin./(out1.totalBin);dist1(isnan(dist1)) = 0; 
% rep 2
out2 = positionDistribution(anal2,bins);
dist2 = out2.eventBin./(out2.totalBin);dist2(isnan(dist2)) = 0; 
% combined 
handle = bar(bins,mean([dist1 dist2],2),'histc');shg
patch('Parent',gca,'Vertices',handle.Vertices,'Faces',handle.Faces,...
    'FaceColor',col,...
    'CData',handle.CData);
set(gca,'Box',true,'FontSize',16,'XLim',lim); 
ylabel(xName)
xlabel('Position relative to stripe center (fraction of AP length)')

%% bcd
figNum = 23;
anal1 = bcdAnal1;
anal2 = bcdAnal2;
col = colr{4};
ax23 = figure(figNum);
% rep 1
out1 = positionDistribution(anal1,bins);
dist1 = out1.eventBin./(out1.totalBin);dist1(isnan(dist1)) = 0; 
% rep 2
out2 = positionDistribution(anal2,bins);
dist2 = out2.eventBin./(out2.totalBin);dist2(isnan(dist2)) = 0; 
% combined 
handle = bar(bins,mean([dist1 dist2],2),'histc');shg
patch('Parent',gca,'Vertices',handle.Vertices,'Faces',handle.Faces,...
    'FaceColor',col,...
    'CData',handle.CData);
set(gca,'Box',true,'FontSize',16,'XLim',lim); 
ylabel(xName)
xlabel('Position relative to stripe center (fraction of AP length)')

%% dst
figNum = 25;
anal1 = dAnal1;
anal2 = dAnal2;
col = colr{5};
ax24 = figure(figNum);
% rep 1
out1 = positionDistribution(anal1,bins);
dist1 = out1.eventBin./(out1.totalBin);dist1(isnan(dist1)) = 0; 
% rep 2
out2 = positionDistribution(anal2,bins);
dist2 = out2.eventBin./(out2.totalBin);dist2(isnan(dist2)) = 0; 
% combined 
handle = bar(bins,mean([dist1 dist2],2),'histc');shg
patch('Parent',gca,'Vertices',handle.Vertices,'Faces',handle.Faces,...
    'FaceColor',col,...
    'CData',handle.CData);
set(gca,'Box',true,'FontSize',16,'XLim',lim); 
ylabel(xName)
xlabel('Position relative to stripe center (fraction of AP length)')
%% save
fp = 'C:\Users\thard\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(ax20,[fp 'spotPostionPDFNeutral'],'svg');
saveas(ax21,[fp 'spotPostionPDFWT'],'svg');
saveas(ax22,[fp 'spotPostionPDFZld'],'svg');
saveas(ax23,[fp 'spotPostionPDFBcd'],'svg');
saveas(ax24,[fp 'spotPostionPDFDst'],'svg');

%% Single exp fits to surv freq; (unused)
% surv Freq, center w/o dropouts
timestep = 15;
fignum = 40;
activationParamsFreqSingle = zeros(5,2);
freqCenterDropout = zeros(5,1);
eventNumberCenterDropout = zeros(5,1);
lowTimeCenterDropout = zeros(5,1);
expLengthMat = zeros(5,2);
meanExpLength = 3.2730e+03; % I got this from averaging all ten exp lengths computed in frequency_dwellFli2
for i = 1:5
    try 
%         out = frequency_dwellFli4(dropout(centerCia1{i}),dropout(centerCia2{i}),dropout(centerCia3{i}),dropout(centerCia4{i}),timestep,colr{i},fignum,[0.1]); %
        out = frequency_dwellFli4((centerCia1{i}),(centerCia2{i}),(centerCia3{i}),(centerCia4{i}),timestep,colr{i},fignum,[0.1]); % 
    catch
%         out = frequency_dwellFli4(dropout(centerCia1{i}),dropout(centerCia2{i}),timestep,colr{i},fignum,[0.1]);
        out = frequency_dwellFli4((centerCia1{i}),(centerCia2{i}),timestep,colr{i},fignum,[0.1]); % 30
    end
    activationParamsFreqSingle(i,:) = [out.totalFrequency out.k1];
end
set(gca,'YLim',[1e-5 1.5e-2],'XLim',[0 2300]);

%% active trxn survival w/ single exp fit; SFig. 6
fignum = 41;

inarg = [0.1];
activationParamsSingle = zeros(5,3);
for i = 1:5
    [out,p] = survivalPlotFromVector(centerCia{i},6,1,1,3600,fignum,inarg,colr{i}); %
    set(p,'Color',colr{i});
    activationParamsSingle(i,:) = [out.tConst out.k1 out.avg];
    hold on;
end
set(gca,'Box',true,'FontSize',16,'YLim',[-6.5 0],'XLim',[0 2300]);
%% save
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
% saveas(gca,[fp 'onTimesFreqAllSingleExpFit220329'],'svg');
% saveas(gca,[fp 'onTimesFreqAllSingleExpFit220401'],'svg'); % with dropouts
saveas(gca,[fp 'onTimesFreqAllSingleExpFit220413'],'svg'); % survival curves

%% SFig 5 raster plots
rasterater4(nAnal1,nAnal2,nAnal3,nAnal4,0,28);
set(gca,'FontSize',16);
rasterater4(spaceAnal1,spaceAnal2,0,29);
set(gca,'FontSize',16);
rasterater4(zldAnal1,zldAnal2,0,30);
set(gca,'FontSize',16);
rasterater4(bcdAnal1,bcdAnal2,0,31);
set(gca,'FontSize',16);
rasterater4(dAnal1,dAnal2,0,32);
set(gca,'FontSize',16);
%% save
figNum = 28;
axN = figure(figNum);axWT = figure(figNum + 1);axZld = figure(figNum + 2);axBcd = figure(figNum + 3);axDst = figure(figNum + 4);
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(axN,[fp 'rasterN2203'],'svg');
saveas(axWT,[fp 'rasterWT2203'],'svg');
saveas(axZld,[fp 'rasterZld2203'],'svg');
saveas(axBcd,[fp 'rasterBcd2203'],'svg');
saveas(axDst,[fp 'rasterDst2203'],'svg');


%% single exp fits to 1st passage SFig 5
% fit time to 1st bind with LJF's scripts (set up):

% see word doc in C:\Users\tth12\matlab\mFiles\timeToFirstBinding
% inputs (from the help menu of exp1_global... .m):
% inarg == [ap r1 r2 r3 ...rN]  where 
%         (active fraction) = ap^2/(1+ap^2)   and
%     r1, r2 r3 are the rates for the rising exponential fits of the form:
%  (active fraction)*(1 - exp(-t*ri) ) - (1-afc)*(1 - exp(-t*rc) )   (2nd term is nonspecific 
%                                               backgrgound surface binding)
% datacells == cell array, where each member contains one experimental data set  
%     i.e. i = 1 to N for N data sets where
% datacells{i}.intervals==vector list of first binding intervals measured for
%                   those AOIS in which at least one landing event was observed
% datacells{i}.tx == maximum observation time 
% datacells{i}.Nt == total number of AOIs, including those already
%                 occupied at t = 0 and those for which no landing was observed
% datacells{i}.Ns == number of AOIs for which no landing was observed
% datacell{i}.Nz == number of AOIs already occupied at time = 0
% rc == rate for nonspecific surface binding 
% call via fminsearch('exp1_global_active_fraction_nonspecific_background_mxl',inarg,[],datacells,rc)

% for i = 1:5
%     %grab the first event (its length)
%     mat = ciaCell{i};
%     first = mat(mat(:,2) == -3,7);
%     %make the first event occur at t = 0;
%     first = first - min(first);
%     dataCells{i}.intervals = first;
%     %mx obs time (you can change this for each exp, but I don't think it
%     %matters):
%     dataCells{i}.tx = 3600;
%     %tot num of nukes. play with this as it has implications to the model
%     dataCells{i}.Nt = NtV(i);
%     %number of nukes with no landing
%     dataCells{i}.Ns = NsV(i)
%     %aois occupied at t = 0:
%     dataCells{i}.Nz = 0; %this is either 0 or 1 if we make first event start at t = 0 
% end
%nonSpecific bind rate. set to zero for the moment. can control for this later
rc = 0;
%% fit time to 1st bind with LJF's scripts (SFig. 5A)
fignum = 67;
NtV = zeros(1,5);
NsV = zeros(1,5);
for i = 1:5
    [out,p] = initialFractionBin(analCell1{i},analCell2{i},bins,boiC,[],30,3000,fignum);
    set(p,'Color',colr{i},'LineWidth',1);
    hold on
    NtV(i) = out.totNukes;
    NsV(i) = NtV(i) - size(out.ints,1); % see below
end

for i = 1:5
    %grab the first event (its length)
    mat = ciaCell{i};
    first = mat(mat(:,2) == -3,7);
    %make the first event occur at t = 0;
    first = first - min(first);
    dataCells{i}.intervals = first;
    %mx obs time (you can change this for each exp, but I don't think it
    %matters):
    dataCells{i}.tx = 3600;
    %tot num of nukes. play with this as it has implications to the model
    dataCells{i}.Nt = NtV(i);
    %number of nukes with no landing
    dataCells{i}.Ns = NsV(i);
    %aois occupied at t = 0:
    dataCells{i}.Nz = 0; %this is either 0 or 1 if we make first event start at t = 0 
end
rc = 0;

t = 1:3600;
inarg = [2 0.001];
kV = [];
afV = [];
for i = 1:5
    outs = fminsearch('exp1_global_active_fraction_nonspecific_background_mxl',inarg,[],dataCells(i),rc); %format of [AF wt zld bcd stat wtSpace];
    afV = [afV; outs(1)^2/(1 + outs(1)^2) ]; %record these for later inspection
    kV = [kV; outs(2)];
    fit = afV(i)*(1 - exp(-outs(2)*t)); %normalized, pair with "A"
%     fit = (afV(i)*dataCells{i}.Nt - dataCells{i}.Nz)*(1 - exp(-(kV(i)+rc)*t)) + dataCells{i}.Nz; %not normalized, pair with "B"
%     ( (Af*Nt - Nz)/Nt )*(rate+rc)*exp(-(rate+rc)*intervals) + (1-Af)*rc*exp(-rc*intervals)
    hold on;plot(t,fit,'--','Color',colr{i});shg
end
set(gca,'YLim',[0 1.01],'Xlim',[0 4000],'Box',true,'FontSize',16); 
ylabel('Fraction of nuclei')
xlabel('Time (s)')

%% SFig 5B: two equal step exp fits to 1st passage
% here we (thought) we did it right. This accounts for limited observation window and
% active fraction of nukes. the rate doesn't change, even with changing
% active fraction
repNum = 5; % for plotting: 1 = neitral, 2 = wt, 3 = zld, 4 = bcd, 5 = dst
% norm = size(spaceCia(spaceCia(:,2) == -3,:),1);  % WT spacer (238 events)
norm = size(zldCia(zldCia(:,2) == -3,:),1);  % Zld (286 events)
fitInputs = [0.1 100]; 
fitMat = [];
% for i = [1 2 4 5] %TSing, WT as norm
for i = 1:repNum
    [fitParams,p,p2] = initialFraction7(centerCia{i},nukeCell{i},NtV(i),fitInputs,30,3000,fignum + 1);
    set(p2,'Color',colr{i});
    fitMat = [fitMat;fitParams.Af fitParams.k0];
%     [out,p] = initialFractionBin(analCell1{i},analCell2{i},bins,boiC,fitInputs,30,3000,fignum + 1);
    set(p,'Color',colr{i},'LineWidth',1);
    hold on;
end
set(gca,'Box',true,'FontSize',16,'XLim',[0 4000],'YLim',[0 1.01]);
ylabel('Fraction of nuclei')
xlabel('Time (s)')
%% save
axOne = figure(fignum);axTwo = figure(fignum + 1);
fp = 'C:\Users\thard\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(axOne,[fp 'firstPassageOneExpFit220328'],'svg');
saveas(axTwo,[fp 'firstPassageTwoEqualExpFit220328'],'svg');










































