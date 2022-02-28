% here is the script to make the Idle distribution SFig
%% load workspace
% neutral spacer eveS2
% 190712 
load t1p51o.dat -mat %this one includes no spot analysis
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
load t1p59e2.dat -mat
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

%% idle time distribution; center, dropout, include time after last event:
%survival, not frequency dist.
% initialize:
fignum = 570;
numBS = 1000;
flag = 0; % 1 to include bootstraps, 0 to not
params = {};
idleParamsCenterDropout = zeros(5,3);
for i = 1:5
    interEventsStruc = idleSurvival2(centerCia{i});
    
    % her are the different forms to plot:
    interEvents = [interEventsStruc.interEvents; interEventsStruc.endEvents]; % use this to include time after last event
    % interEvents = [interEventsStruc.interEvents]; % no last events.
    
    % fits:
    % [params,p] = survivalPlotFromVector(interEvents,1,0,1,3600,fignum); % no fit
    [params{i},p,f] = survivalPlotFromVector2(interEvents,1,2,1,3600,fignum,[0.8 0.01 0.001]); % dbl exp
%     [params{i},p] = survivalPlotFromVector(interEvents,1,1,1,3600,fignum,0.1); % single exp
    % b = params.a;u1 = params.k1;u2 = params.k2;
    set(p,'Color',colr{i});
    set(f,'Color',colr{i},'LineStyle','--','LineWidth',1);
    hold on
    if flag == 1
        bootstrapsBootstraps(interEvents,1,numBS,colr{i},fignum + 1);
        hold on
    end
    % store the fit values:
    idleParamsCenterDropout(i,:) = [params{i}.a params{i}.k1 params{i}.k2];
end
set(gca,'Box',true,'FontSize',16,'YLim',[-6 0]);
% set(gca,'YScale','log');
% ylabel('Ln(probability of survival)')
% xlabel('Time (s)')
% off times doesn't work, as in it is unfittable, if we include the times
% after the last event. So I think this may be useless. ...220215: this is
% wrong

% 220215 update:
% these are are very fitable.
% Dbl exp Fits of the inter event times only yield something close to
% single exponential results (all the As are close to one.

%% idle time bootstraps; center, dropout, include time after last event:
inarg = [1 3600 0.8 0.01 0.001]; %[tm tx B u_1 u_2]
numBS = 1000;
constructs = 1:5;
idleCenterFits = {};
idleParamsCenterDropoutStd1000 = zeros(5,3);
tic;
for i = constructs
    % get the idle events (again)
    interEventsStruc = idleSurvival2(centerCia{i});
    interEvents = [interEventsStruc.interEvents; interEventsStruc.endEvents];
    % now boot strap and fit:
    bsFits = bootstrapsBootstrapsFit(interEvents,1,numBS,2,inarg);
    idleCenterFits{i} = bsFits;
    idleParamsCenterDropoutStd1000(i,:) = bsFits.std;
end
toc

%% table it (including time after last event):
% set up table
Reporter = tits';
% rates:
% allocate:
A = zeros(5,1);k1 = A; k2 = A; %allocate space, amke things col vectors
Aerr = A; k1err = A; k2err = A;
errMat = idleParamsCenterDropoutStd1000;
for i = 1:5
    A(i) = round(idleParamsCenterDropout(i,1),2);
    Aerr(i) = round(errMat(i,1),2);
    k1(i) = (idleParamsCenterDropout(i,2));
    k1err(i) = (errMat(i,2));
    k2(i) = (idleParamsCenterDropout(i,3));
    k2err(i) = (errMat(i,3));
end
% taus:
t1 = 1./k1;
t2 = 1./k2;
%propagate errror:
t1err = round(1./k1.^2.*k1err,-1);
t2err = round(1./k2.^2.*k2err,-1);
% round
t1 = round(t1,-1);
t2 = round(t2,-1);
activeParamsCenterDropout = table(Reporter,A,Aerr,t1,t1err,t2,t2err)

%% idle time distribution; center, dropout, NO time after last event:
%initialize:
fignum = 72;
numBS = 1000;
flag = 0; % 1 to include bootstraps, 0 to not
params = {};
idleParamsCenterDropout2 = zeros(5,3);
for i = 1:5
    interEventsStruc = idleSurvival2(centerCia{i});
    
    % her are the different forms to plot:
    interEvents = [interEventsStruc.interEvents]; % use this to include time after last event
    % interEvents = [interEventsStruc.interEvents]; % no last events.
    
    % fits:
    % [params,p] = survivalPlotFromVector(interEvents,1,0,1,3600,fignum); % no fit
    [params{i},p,f] = survivalPlotFromVector2(interEvents,1,2,1,3600,fignum,[0.8 0.01 0.001]); % dbl exp
%     [params{i},p] = survivalPlotFromVector(interEvents,1,1,1,3600,fignum,0.1); % single exp
    % b = params.a;u1 = params.k1;u2 = params.k2;
    set(p,'Color',colr{i});
    set(f,'Color',colr{i},'LineStyle','--','LineWidth',1);
    hold on
    if flag == 1
        bootstrapsBootstraps(interEvents,1,numBS,colr{i},fignum);
    end
    % store the fit values:
    idleParamsCenterDropout2(i,:) = [params{i}.a params{i}.k1 params{i}.k2];
end
set(gca,'Box',true,'FontSize',16,'YLim',[-6 0]);

%% idle time bootstraps; center, dropout, NO time after last event:
inarg = [1 3600 0.8 0.01 0.001]; %[tm tx B u_1 u_2]
constructs = 1:5;
idleCenterFits2 = {};
idleParamsCenterDropout2Std1000 = zeros(5,3);
tic;
for i = constructs
    % get the idle events (again)
    interEventsStruc = idleSurvival2(centerCia{i});
    interEvents = [interEventsStruc.interEvents]; % no end events
    % now boot strap and fit:
    bsFits = bootstrapsBootstrapsFit(interEvents,1,numBS,2,inarg);
    idleCenterFits2{i} = bsFits;
    idleParamsCenterDropout2Std1000(i,:) = bsFits.std;
end
toc

%% table it (NO time after last event):
% set up table
Reporter = tits';
% rates:
% allocate:
A = zeros(5,1);k1 = A; k2 = A; %allocate space, amke things col vectors
Aerr = A; k1err = A; k2err = A;
errMat = idleParamsCenterDropout2Std1000;
for i = 1:5
    A(i) = round(idleParamsCenterDropout2(i,1),2);
    Aerr(i) = round(errMat(i,1),2);
    k1(i) = (idleParamsCenterDropout2(i,2));
    k1err(i) = (errMat(i,2));
    k2(i) = (idleParamsCenterDropout2(i,3));
    k2err(i) = (errMat(i,3));
end
% taus:
t1 = 1./k1;
t2 = 1./k2;
%propagate errror:
t1err = round(1./k1.^2.*k1err,-1);
t2err = round(1./k2.^2.*k2err,-1);
% round
t1 = round(t1,-1);
t2 = round(t2,-1);
activeParamsCenterDropout = table(Reporter,A,Aerr,t1,t1err,t2,t2err)

%% idle time BAR CHART; center, dropout, NO time after last event:
% first calc the frequencies:

timestep = 15;
fignum = 231;
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
        out = frequency_dwellFli5(dropout(centerCia1{i}),dropout(centerCia2{i}),dropout(centerCia3{i}),dropout(centerCia4{i}),timestep,colr{i},fignum,[0.9 0.01 0.001],[]); % 
            expLengthMat(i,:) = out.expLength;
    catch
        out = frequency_dwellFli5(dropout(centerCia1{i}),dropout(centerCia2{i}),timestep,colr{i},fignum,[0.9 0.01 0.001],[]); % 
        expLengthMat(i,1:2) = out.expLength;
    end
    activationParamsCenterDropout(i,:) = [out.a out.k1 out.k2];
    freqCenterDropout(i) = 1/out.totalFrequency;
    eventNumberCenterDropout(i) = out.eventNumber;
    lowTimeCenterDropout(i) = out.lowTimes;
end
set(gca,'YLim',[7e-6 5e-3],'XLim',[0 2300]);

% calc the idle periods, then spit them out in a table:
Reporter = tits';
% counting error:
eventErr = sqrt(eventNumberCenterDropout); 
% frequencies:
f = round(eventNumberCenterDropout./lowTimeCenterDropout,4);
% error:
ferr = round(eventErr./lowTimeCenterDropout,4);
% period: 
T = round(freqCenterDropout,-1);
% error
% Terr = round(lowTimeCenterDropout./eventNumberCenterDropout.^2.*eventErr,-1);
Terr = round(eventErr,-1);
freqCenterDropoutTable = table(Reporter,f,ferr,T,Terr)
%% then make the bar chart:
figure(230);
bar(1,T(1),0.4,'FaceColor',colr{1},'EdgeColor','none'); %initiate the bar chart first
for i = 2:5
    hold on; bar(i,T(i),0.4,'FaceColor',colr{i},'EdgeColor','none');
end
xticks([1 2 3 4 5]);
set(gca,'xticklabels',tits(1:5),'yLim',[0 450],'Box',true,'FontSize',16,'xLim',[0.5 5.5]);
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



































% idleSurvivalNoDropoutCenterDblFitErrEndEvents220218