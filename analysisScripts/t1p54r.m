% to produce the image acquiosition conditions figure (i.e. SFig 1)
% of Harden et al 2022

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

Reporter = tits';

%% surv Freq, center, all times, (need for fig. 3E)
timestep = 15;
fignum = 30;
activationParamsCenterDropout = zeros(5,3);
freqCenterDropout = zeros(5,1);
eventNumberCenterDropout = zeros(5,1);
lowTimeCenterDropout = zeros(5,1);
expLengthMat = zeros(5,2);
meanExpLength = 3.2730e+03; % I got this from averaging all ten exp lengths computed in frequency_dwellFli2
for i = 1:5
    try
        out = frequency_dwellFli5(centerCia1{i},centerCia2{i},centerCia3{i},centerCia4{i},timestep,colr{i},fignum,[0.9 0.01 0.001],[]); % 30
        activationParamsCenterDropout(i,:) = [out.a out.k1 out.k2];
        expLengthMat(i,:) = out.expLength;
        freqCenterDropout(i) = 1/out.totalFrequency;
        eventNumberCenterDropout(i) = out.eventNumber;
        lowTimeCenterDropout(i) = out.lowTimes;
    catch
        out = frequency_dwellFli5(centerCia1{i},centerCia2{i},timestep,colr{i},fignum,[0.9 0.01 0.001],[]); % 30
        activationParamsCenterDropout(i,:) = [out.a out.k1 out.k2];
        expLengthMat(i,:) = out.expLength;
        freqCenterDropout(i) = 1/out.totalFrequency;
        eventNumberCenterDropout(i) = out.eventNumber;
        lowTimeCenterDropout(i) = out.lowTimes;
    end
end
set(gca,'YLim',[1e-5 1.5e-2],'XLim',[0 2300]);

%% idle period, table, all times, Table 3
% define some terms: 
% 1. event count: eventNumberCenterDropout
% 2. total low time: lowTimeCenterDropout
Reporter = tits';
% counting error:
eventErr = sqrt(eventNumberCenterDropout); 
% frequencies:
f = round(eventNumberCenterDropout./lowTimeCenterDropout,4);
% error:
ferr = round(eventErr./lowTimeCenterDropout,4);
% period: 
T = round(freqCenterDropout);
% error
Terr = round(lowTimeCenterDropout./eventNumberCenterDropout.^2.*eventErr);
freqCenterDropoutTable = table(Reporter,f,ferr,T,Terr)
%% idle period bar chart, all times Fig 3E
ax220 = figure(220);
bar(1,T(1),0.4,'FaceColor',colr{1},'EdgeColor','none'); %initiate the bar chart first
for i = 2:5
    hold on; bar(i,T(i),0.4,'FaceColor',colr{i},'EdgeColor','none');
end
xticks([1 2 3 4 5]);
set(gca,'xticklabels',tits(1:5),'yLim',[0 200],'Box',true,'FontSize',16,'xLim',[0.5 5.5]);
xticklabel_rotate([],45,[],'FontSize',16);
hold on;
er = errorbar(1,T(1),Terr(1),Terr(1));
er.Color = colr{2};    %make the black neutral grey                       
er.LineStyle = 'none'; 
for i = 2:5
    hold on;
    er = errorbar(i,T(i),Terr(i),Terr(i));
    er.Color = colr{1};    % make all err bars black                       
    er.LineStyle = 'none'; 
end
set(gca,'Box',true);
%% save 
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(gca,[fp 'idlePeriodAllTimes220506'],'svg');

%% surv Freq, EARLY, center, unused for figures
%used in paper figs
fignum = 34;
split = 1/5; %Sets the divide btwn the two bins set to 1/2 for 50/50, 1/3 for first third in "early" bin, etc
timestep = 15;
onHalf1cell1 = {};onHalf1cell2 = {};freqCenterDropoutEarly = [];
activationParamsCenterDropoutEarly = zeros(5,3);
freqCenterDropoutEarly = zeros(5,1);
eventNumberCenterDropoutEarly = zeros(5,1);
lowTimeCenterDropoutEarly = zeros(5,1);
flag = 0; % 1 for bootstraps, 0 for not
constructs = 1:5;
for i = constructs
    % treat each rep indie:
    mat1 = (centerCia1{i});
    L = size(mat1,1);
    ts = sort(mat1(:,7));
    t50 = ts(floor(L*split));
    logi = mat1(:,7) <= t50;
    half1mat = mat1(logi,:);
    onHalf1cell1{i} = (half1mat);
    % the other rep:
    mat2 = (centerCia2{i});
    L = size(mat2,1);
    ts = sort(mat2(:,7));
    t50 = ts(floor(L*split));
    logi = mat2(:,7) <= t50;
    half1mat = mat2(logi,:);
    onHalf1cell2{i} = (half1mat);
    try
        % 3rd rep:
        mat3 = (centerCia3{i});
        L = size(mat3,1);
        ts = sort(mat3(:,7));
        t50 = ts(floor(L*split));
        logi = mat3(:,7) <= t50;
        half1mat = mat3(logi,:);
        onHalf1cell3{i} = (half1mat);
        % 4th rep:
        mat4 = (centerCia4{i});
        L = size(mat4,1);
        ts = sort(mat4(:,7));
        t50 = ts(floor(L*split));
        logi = mat4(:,7) <= t50;
        half1mat = mat4(logi,:);
        onHalf1cell4{i} = (half1mat);
        [out,p] = frequency_dwellFli5(onHalf1cell1{i},onHalf1cell2{i},onHalf1cell3{i},onHalf1cell4{i},timestep,colr{i},fignum,[0.9 0.01 0.001],[]);
%         [out,p] = survivalPlotFromVector([onHalf1cell1{i};onHalf1cell2{i};onHalf1cell3{i};onHalf1cell4{i}],6,2,1,3600,fignum,[0.9 0.01 0.001],colr{i});
%         set(p,'Color',colr{i});
        activationParamsCenterDropoutEarly(i,:) = [out.a out.k1 out.k2];
        freqCenterDropoutEarly(i) = 1/out.totalFrequency;
        eventNumberCenterDropoutEarly(i) = out.eventNumber;
        lowTimeCenterDropoutEarly(i) = out.lowTimes;
    catch
        [out,p] = frequency_dwellFli5(onHalf1cell1{i},onHalf1cell2{i},timestep,colr{i},fignum,[0.9 0.01 0.001],[]);
%         [out,p] = survivalPlotFromVector([onHalf1cell1{i};onHalf1cell2{i}],6,2,1,3600,fignum,[0.9 0.01 0.001],colr{i}); %
%         set(p,'Color',colr{i});
        activationParamsCenterDropoutEarly(i,:) = [out.a out.k1 out.k2];
        freqCenterDropoutEarly(i) = 1/out.totalFrequency;
        eventNumberCenterDropoutEarly(i) = out.eventNumber;
        lowTimeCenterDropoutEarly(i) = out.lowTimes;
        hold on
    end
    if flag == 1
        try
            [~,pbs] = bootstrapsFrequency([onHalf1cell1{i};onHalf1cell2{i};onHalf1cell3{i};onHalf1cell4{i}],6,10000,freqCenterDropoutEarly(i),colr{i},fignum);
            hold on
        catch 
            [~,pbs] = bootstrapsFrequency([onHalf1cell1{i};onHalf1cell2{i}],6,10000,freqCenterDropoutEarly(i),colr{i},fignum);
            hold on
        end
    end
end
%plot stuff:
% set(gca,'YLim',[2e-4 2.5e-2],'XLim',[0 2400]); % neutral & wt
set(gca,'YLim',[2e-4 0.06],'XLim',[0 2400]); % all cndns

%% idel period table,  EARLY, Table 3
eventErrEarly = sqrt(eventNumberCenterDropoutEarly); 
earlyf = round(eventNumberCenterDropoutEarly./lowTimeCenterDropoutEarly,3);
earlyferr = round(eventErrEarly./lowTimeCenterDropoutEarly,3);
earlyT = round(freqCenterDropoutEarly);
earlyTerr = round(lowTimeCenterDropoutEarly./eventNumberCenterDropoutEarly.^2.*eventErrEarly);
freqCenterDropoutTableEarly = table(Reporter,earlyf,earlyferr,earlyT,earlyTerr)
%% idel period bar EARLY
ax221 = figure(221);
bar(1,earlyT(1),0.4,'FaceColor',colr{1},'EdgeColor','none'); %initiate the bar chart first
for i = 2:5
    hold on; bar(i,earlyT(i),0.4,'FaceColor',colr{i},'EdgeColor','none');
end
xticks([1 2 3 4 5]);
set(gca,'xticklabels',tits(1:5),'yLim',[0 55],'Box',true,'FontSize',16,'xLim',[0.5 5.5]);
xticklabel_rotate([],45,[],'FontSize',16);
hold on;
er = errorbar(1,earlyT(1),earlyTerr(1),earlyTerr(1));
er.Color = colr{2};    %make the black neutral grey                       
er.LineStyle = 'none'; 
for i = 2:5
    hold on;
    er = errorbar(i,earlyT(i),earlyTerr(i),earlyTerr(i));
    er.Color = colr{1};    % make all err bars black                       
    er.LineStyle = 'none'; 
end
set(gca,'Box',true);

%% save 
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(gca,[fp 'idlePeriodEarly220506'],'svg');

%% surv Freq, LATE, center, unused for figures
% not the idle periods between the late and early cias are getting lost,
% let's recover them and stick them into the late distributions:
% must run this for the EARLY events first
fignum = 35;
timestep = 15;
onHalf2cell1 = {};onHalf2cell2 = {};onLateParamsCell = {};
activationParamsCenterDropoutLate = zeros(5,3);
freqCenterDropoutLate = zeros(5,1);
eventNumberCenterDropoutLate = zeros(5,1);
lowTimeCenterDropoutLate = zeros(5,1);
flag = 0; 
for i = 1:5
    % treat each rep indie:
    mat1 = centerCia1{i};
    u = unique(mat1(:,1));
    half2mat = [];
    for j = 1:size(u,1)
        nuke = u(j);
        % get all the events for this nuke:
        nukeMat = mat1(mat1(:,1) == nuke,:);
        % how many events are there?
        eventsInThisNuke = size(nukeMat,1);
        % check to see if this nuke has early events:
        logi = onHalf1cell1{i}(:,1) == nuke;
        if sum(logi) > 0
            % get the cia for early events in this nuke:
            earlyNukeMat = onHalf1cell1{i}(logi,:);
            % how many early events are in this nuke?
            earlyEventsInThisNuke = size(earlyNukeMat,1);
            % if this nuke has early and late events:
            if eventsInThisNuke > earlyEventsInThisNuke
                % get the last event included in the early data set and all
                % other events:
                addOn = nukeMat(earlyEventsInThisNuke:end,:);
                half2mat = [half2mat; addOn];
            end
        % if this nuke doesn't ahve early events, all its events must be
        % late:
        else
            half2mat = [half2mat; nukeMat];
        end
    end
    onHalf2cell1{i} = half2mat;
    % the other rep:
    mat2 = centerCia2{i};
    u = unique(mat2(:,1));
    half2mat = [];
    for j = 1:size(u,1)
        nuke = u(j);
        nukeMat = mat2(mat2(:,1) == nuke,:);
        eventsInThisNuke = size(nukeMat,1);
        logi = onHalf1cell2{i}(:,1) == nuke;
        if sum(logi) > 0
            earlyNukeMat = onHalf1cell2{i}(logi,:);
            earlyEventsInThisNuke = size(earlyNukeMat,1);
            if eventsInThisNuke > earlyEventsInThisNuke
                addOn = nukeMat(earlyEventsInThisNuke:end,:);
                half2mat = [half2mat; addOn];
            end
        else
            half2mat = [half2mat; nukeMat];
        end
    end
    onHalf2cell2{i} = half2mat;
    try
        % 3rd rep:
        mat3 = centerCia3{i};
        u = unique(mat3(:,1));
        half2mat = [];
        for j = 1:size(u,1)
            nuke = u(j);
            nukeMat = mat3(mat3(:,1) == nuke,:);
            eventsInThisNuke = size(nukeMat,1);
            logi = onHalf1cell3{i}(:,1) == nuke;
            if sum(logi) > 0
                earlyNukeMat = onHalf1cell3{i}(logi,:);
                earlyEventsInThisNuke = size(earlyNukeMat,1);
                if eventsInThisNuke > earlyEventsInThisNuke
                    addOn = nukeMat(earlyEventsInThisNuke:end,:);
                    half2mat = [half2mat; addOn];
                end
            else
                half2mat = [half2mat; nukeMat];
            end
        end
        onHalf2cell3{i} = (half2mat);
        % 4th rep:
        mat4 = centerCia4{i};
        u = unique(mat4(:,1));
        half2mat = [];
        for j = 1:size(u,1)
            nuke = u(j);
            nukeMat = mat4(mat4(:,1) == nuke,:);
            eventsInThisNuke = size(nukeMat,1);
            logi = onHalf1cell4{i}(:,1) == nuke;
            if sum(logi) > 0
                earlyNukeMat = onHalf1cell4{i}(logi,:);
                earlyEventsInThisNuke = size(earlyNukeMat,1);
                if eventsInThisNuke > earlyEventsInThisNuke
                    addOn = nukeMat(earlyEventsInThisNuke:end,:);
                    half2mat = [half2mat; addOn];
                end
            else
                half2mat = [half2mat; nukeMat];
            end
        end
        onHalf2cell4{i} = half2mat;
        [out,p] = frequency_dwellFli5(onHalf2cell1{i},onHalf2cell2{i},onHalf2cell3{i},onHalf2cell4{i},timestep,colr{i},fignum,[0.9 0.01 0.001],[]);
        activationParamsCenterDropoutLate(i,:) = [out.a out.k1 out.k2];
        freqCenterDropoutLate(i) = 1/out.totalFrequency;
        eventNumberCenterDropoutLate(i) = out.eventNumber;
        lowTimeCenterDropoutLate(i) = out.lowTimes;
        hold on
        set(p,'Color',colr{i});
    catch
        [out,p] = frequency_dwellFli5(onHalf2cell1{i},onHalf2cell2{i},timestep,colr{i},fignum,[0.9 0.01 0.001],[]);
        activationParamsCenterDropoutLate(i,:) = [out.a out.k1 out.k2];
        freqCenterDropoutLate(i) = 1/out.totalFrequency;
        eventNumberCenterDropoutLate(i) = out.eventNumber;
        lowTimeCenterDropoutLate(i) = out.lowTimes;
        hold on
        set(p,'Color',colr{i});
    end
    if flag == 1
        [~,pbs] = bootstrapsFrequency([onHalf2cell1{i};onHalf2cell2{i}],6,1000,freqCenterDropoutLate(i),colr{i},fignum);
        hold on
    end
end
set(gca,'YLim',[1.0e-5 0.006],'XLim',[0 1700]); % all cndns

%% idel period, table LATE, Table 3
eventErrLate = sqrt(eventNumberCenterDropoutLate); 
latef = round(eventNumberCenterDropoutLate./lowTimeCenterDropoutLate,4);
lateferr = round(eventErrLate./lowTimeCenterDropoutLate,4);
lateT = round(freqCenterDropoutLate);
lateTerr = round(lowTimeCenterDropoutLate./eventNumberCenterDropoutLate.^2.*eventErrLate);
freqCenterDropoutTableLate = table(Reporter,latef,lateferr,lateT,lateTerr)

%% idel period bar LATE, unused for figures
ax221 = figure(222);
bar(1,lateT(1),0.4,'FaceColor',colr{1},'EdgeColor','none'); %initiate the bar chart first
for i = 2:5
    hold on; bar(i,lateT(i),0.4,'FaceColor',colr{i},'EdgeColor','none');
end
xticks([1 2 3 4 5]);
set(gca,'xticklabels',tits(1:5),'yLim',[0 250],'Box',true,'FontSize',16,'xLim',[0.5 5.5]);
xticklabel_rotate([],45,[],'FontSize',16);
hold on;
er = errorbar(1,lateT(1),lateTerr(1),lateTerr(1));
er.Color = colr{2};    %make the black neutral grey                       
er.LineStyle = 'none'; 
for i = 2:5
    hold on;
    er = errorbar(i,lateT(i),lateTerr(i),lateTerr(i));
    er.Color = colr{1};    % make all err bars black                       
    er.LineStyle = 'none'; 
end
set(gca,'Box',true);
