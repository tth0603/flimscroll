% to produce the image acquiosition conditions figure (i.e. SFig 1)
% of Harden et al 2022

%% initialize variables:
% 1/30s exposure (2 Frame/min)
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

% 1/60s (1 frsme/min)
% 211027
load t1p71g2.dat -mat
space60anal = analData;
clear analData
space60cia = space60anal.cia;
space60spots = space60anal.spotsMat;
space60nukes = space60anal.nukesMat;

% 1/15s (4 frame/min)
% 211026
load t1p71l2.dat -mat
space15anal = analData;
clear analData
space15cia = space15anal.cia;
space15spots = space15anal.spotsMat;
space15nukes = space15anal.nukesMat;

% line color set up
clear colr
colr{1} = [0 0 0]; %15s blk (4 frame/min)
colr{2} = [0.6 0 0]; %30s dark red(2 frame/min)
colr{3} = [0.6 0.6 0.6];  %60s dark grey (1 frame/min)

% legend set up 
clear tits
tits{1} = 'eve2:wt; 4 frame/min';
tits{2} = 'eve2:wt; 2 frame/min';
tits{3} = 'eve2:wt; 1 frame/min';

% nuke binning stuff
medNukeNum = 154.5; % this was empirically determined in t1p54l.m, so we'll use it here for consistency.
% this only matters for spots v T, which we likely won't use

% cells for loops:
% analysis cell:
analCell{1} = space15anal;
analCell{2} = spaceAnal;
analCell{3} = space60anal;
% cic cell:
ciaCell{1} = space15cia;
ciaCell{2} = spaceCia;
ciaCell{3} = space60cia;

% plotting stuff:
dutyCycle = [4 2 1]'; % [15s 30s 60s] conditions (in frame/min)

% make cells limit analysis to strip center:
bins = [-0.20:0.02:0.20]; % define the bins
boiC = 10:11; % for the center bins

%take care of neutral first bc of the extra reps
centerCia{1} = binnedCiaMaker3(space15anal,bins,boiC);
centerCia{2} = binnedCiaMaker3(spaceAnal1,spaceAnal2,bins,boiC);
centerCia1{1} = binnedCiaMaker3(spaceAnal1,bins,boiC);
centerCia2{1} = binnedCiaMaker3(spaceAnal2,bins,boiC);
centerCia{3} = binnedCiaMaker3(space60anal,bins,boiC);

% idle cias:
idleCia = {};idleCia1 = {};idleCia2 = {};
for i = 1:3
    idleCia{i} = idleCiaMaker2(centerCia{i},30);
    try
        idleCia1{i} = idleCiaMaker2(centerCia1{i},30);
        idleCia2{i} = idleCiaMaker2(centerCia2{i},30);
    end
end

%% spots v T (SFig. 1A)
% 2 frame/min:
figure(535);
clear cia15dropout mat2 mn mx
%get all rose with a spot:
cia15dropout = spaceNukes(spaceNukes(:,4) > 0,:); %col 4 is binary
%get frame range:
mx = max(cia15dropout(:,2)); %frame indexing is more useful than time
svt = []; %spots v time
for i = 1:mx
    mat2 = cia15dropout(cia15dropout(:,2) == i,:); %find all rose for this frame (logical indexing)
    if ~isempty(mat2) %if we have nukes with spots, mark the time and sum up the num of spots
        svt = [svt;mat2(1,1) size(mat2,1)]; %[time numOfSpots]
    end
end
hold on;plot(svt(:,1),svt(:,2)./medNukeNum,'Color',colr{2},'LineStyle','-','MarkerFaceColor',colr{2},'MarkerEdgeColor','none','MarkerSize',3,'Marker','o');shg

% 1 frame/min:
clear cia15dropout mat2 mn mx
%get all rose with a spot:
cia15dropout = space60nukes(space60nukes(:,4) > 0,:); %col 4 is binary
%get frame range:
mx = max(cia15dropout(:,2)); %frame indexing is more useful than time
svt = []; %spots v time
for i = 1:mx
    mat2 = cia15dropout(cia15dropout(:,2) == i,:); %find all rose for this frame (logical indexing)
    if ~isempty(mat2) %if we have nukes with spots, mark the time and sum up the num of spots
        svt = [svt;mat2(1,1) size(mat2,1)]; %[time numOfSpots]
    end
end
hold on;plot(svt(:,1),svt(:,2)./(medNukeNum/2),'Color',colr{3},'LineStyle','-','MarkerFaceColor',colr{3},'MarkerEdgeColor','none','MarkerSize',3,'Marker','o');shg
% use '(medNukeNum/2)' here because we only have a single replicates, so we
% take the medium number of nukes per bin and divide by two

% 4 frame/min:
clear cia15dropout mat2 mn mx
%get all rose with a spot:
cia15dropout = space15nukes(space15nukes(:,4) > 0,:); %col 4 is binary
%get frame range:
mx = max(cia15dropout(:,2)); %frame indexing is more useful than time
svt = []; %spots v time
for i = 1:mx
    mat2 = cia15dropout(cia15dropout(:,2) == i,:); %find all rose for this frame (logical indexing)
    if ~isempty(mat2) %if we have nukes with spots, mark the time and sum up the num of spots
        svt = [svt;mat2(1,1) size(mat2,1)]; %[time numOfSpots]
    end
end
hold on;plot(svt(:,1),svt(:,2)./(medNukeNum/2),'Color',colr{1},'LineStyle','-','MarkerFaceColor',colr{1},'MarkerEdgeColor','none','MarkerSize',3,'Marker','o');shg
% use '(medNukeNum/2)' here because we only have a single replicates, so we
% take the medium number of nukes per bin and divide by two

set(gca,'Box',true,'FontSize',16,'xlim',[0 4500]);
% ylabel('Number of spots');xlabel('Time (s)');

% NB: This plot looks as if there is no discernible PBing in these
% conditions (a miracle), but there is a bump at about 30min (frame 2000)
% in number of spots in the 1 frame/min condition. 

%% active trxn survival, center, SFig. 1C
fignum = 530;
numOfBS = 10000;
flag = 1;
inarg = [0.9 0.01 0.001];
activationParamsCenterPB = zeros(3,3);
for i = 1:3
    [out,p] = survivalPlotFromVector((centerCia{i}),6,2,1,3600,fignum,inarg,colr{i}); %
    set(p,'Color',colr{i});
    if flag == 1 
        [~,pbs] = bootstrapsBootstraps((centerCia{i}),6,numOfBS,colr{i},fignum);
    end
    hold on;
    activationParamsCenterPB(i,:) = [out.a out.k1 out.k2];
end
set(gca,'Box',true,'FontSize',16,'YLim',[-6.5 0],'XLim',[0 2300]);

%% save this plot:
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(gca,[fp 'survCurveFit&errPB220413'],'svg'); % 

%% surv Freq of active transcription times. CENTER spots. (need for SFig. 1D)
% use this to be consistent
% initialize:
fignum = 531;
activationParamsCenterDropout = zeros(3,3);
freqCenterDropout = zeros(3,1);
eventNumberCenterDropout = zeros(3,1);
lowTimeCenterDropout = zeros(3,1);
expLengthMat = zeros(3,4);
meanExpLength = 3.2730e+03; % I got this from averaging all ten exp lengths computed in frequency_dwellFli2

% 4 frame/min:
timestep = 15; %using half of the acqn time
% out = frequency_dwellFli3(dropout(centerCia{1}),timestep,colr{1},fignum,[0.9 0.01 0.001],[]); % 
out = frequency_dwellFli5((centerCia{1}),timestep,colr{1},fignum,[0.9 0.01 0.001],[]); % 
expLengthMat(1,1:2) = out.expLength;
activationParamsCenterDropout(1,:) = [out.a out.k1 out.k2];
freqCenterDropout(1) = 1/out.totalFrequency;
eventNumberCenterDropout(1) = out.eventNumber;
lowTimeCenterDropout(1) = out.lowTimes;

% 2 frame/min:
timestep = 30; %using half of the acqn time
hold on;
% out = frequency_dwellFli3(dropout(centerCia{2}),timestep,colr{2},fignum,[0.9 0.01 0.001],[]); 
% out = frequency_dwellFli3(dropout(centerCia1{1}),dropout(centerCia2{1}),timestep,colr{2},fignum,[0.9 0.01 0.001],[]); % 
out = frequency_dwellFli5((centerCia1{1}),(centerCia2{1}),timestep,colr{2},fignum,[0.9 0.01 0.001],[]); % % treat the two replicates separately
expLengthMat(2,1:2) = out.expLength;
activationParamsCenterDropout(2,:) = [out.a out.k1 out.k2];
freqCenterDropout(2) = 1/out.totalFrequency;
eventNumberCenterDropout(2) = out.eventNumber;
lowTimeCenterDropout(2) = out.lowTimes;

% 1 frame/min:
timestep = 60; %using half of the acqn time
hold on;
% out = frequency_dwellFli3(dropout(centerCia{3}),timestep,colr{3},fignum,[0.9 0.01 0.001],[]); % 
out = frequency_dwellFli5((centerCia{3}),timestep,colr{3},fignum,[0.9 0.01 0.001],[]); % 
expLengthMat(3,1:2) = out.expLength;
activationParamsCenterDropout(3,:) = [out.a out.k1 out.k2];
freqCenterDropout(3) = 1/out.totalFrequency;
eventNumberCenterDropout(3) = out.eventNumber;
lowTimeCenterDropout(3) = out.lowTimes;
set(gca,'YLim',[9e-6 1e-2],'XLim',[0 2300]);

%% idle period, table, all times, (need for SFig 1D)
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
%% idle period bar chart SFig 1D
ax920 = figure(920);
bar(1,T(1),0.4,'FaceColor',colr{1},'EdgeColor','none'); %initiate the bar chart first
for i = 2:3
    hold on; bar(i,T(i),0.4,'FaceColor',colr{i},'EdgeColor','none');
end
xticks([1 2 3]);
set(gca,'xticklabels',tits(1:3),'yLim',[0 230],'Box',true,'FontSize',16,'xLim',[0.5 3.5]);
xticklabel_rotate([],45,[],'FontSize',16);
hold on;
er = errorbar(1,T(1),Terr(1),Terr(1));
er.Color = colr{3};    %make the black neutral grey                       
er.LineStyle = 'none'; 
for i = 2:3
    hold on;
    er = errorbar(i,T(i),Terr(i),Terr(i));
    er.Color = colr{1};    % make all err bars black                       
    er.LineStyle = 'none'; 
end
set(gca,'Box',true);
%% save 
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(gca,[fp 'idlePeriodPB220516'],'svg');

%% active trxn fit errors (SFig 1C)
inarg = [1 3600 0.9 0.01 0.001]; %[tm tx a k1 k2]
numBs = 10000;
pbActiveTrxnFits = {};
pbActiveTrxnFitsStd = zeros(3,3);
tic;
% 15 s
bsFits = bootstrapsBootstrapsFit((centerCia{1}),6,numBs,2,inarg);
pbActiveTrxnFits{1} = bsFits;
pbActiveTrxnFitsStd(1,:) = bsFits.std;
% 30 s
bsFits = bootstrapsBootstrapsFit((centerCia{2}),6,numBs,2,inarg);
pbActiveTrxnFits{2} = bsFits;
pbActiveTrxnFitsStd(2,:) = bsFits.std;
% 60 s
bsFits = bootstrapsBootstrapsFit((centerCia{3}),6,numBs,2,inarg);
pbActiveTrxnFits{3} = bsFits;
pbActiveTrxnFitsStd(3,:) = bsFits.std;
toc
% should take about 8 min on your lab machine with bs = 10000

%% table. active trxn PBing params (SFig. 1C)
% set up table
Reporter = tits';
% rates:
% allocate:
A = zeros(3,1);k1 = A; k2 = A; %allocate space, amke things col vectors
Aerr = A; k1err = A; k2err = A;
errMat = pbActiveTrxnFitsStd;
for i = 1:3
    A(i) = round(activationParamsCenterDropout(i,1),2);
    Aerr(i) = round(errMat(i,1),2);
    k1(i) = (activationParamsCenterDropout(i,2));
    k1err(i) = round(errMat(i,2),3);
    k2(i) = (activationParamsCenterDropout(i,3));
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
activeParamsDropout = table(Reporter,A,Aerr,t1,t1err,t2,t2err)


% plot the fit params (not used to make figures):
%% amplitude, a; 
t = [0:0.01:5];
figure(506);
errorbar(dutyCycle,A,Aerr,'k.','MarkerSize',15);shg
set(gca,'Ylim',[0 1],'Xlim',[0 5],'Box',true,'FontSize',16);
% add the (weighted) fit:
ampFit = weightedfit([dutyCycle A Aerr]);
% make the line:
yInt = ampFit.Intercept;
m = ampFit.slope;
ampFitCurve = m.*t + yInt;
hold on;plot(t,ampFitCurve,'r');shg
%% tau1;
figure(507);
errorbar(dutyCycle,t1,t1err,'.','MarkerSize',15);shg
set(gca,'Ylim',[0 150],'Xlim',[0 5],'Box',true,'FontSize',16);

%% tau2;
figure(508);
errorbar(dutyCycle,t2,t2err,'.','MarkerSize',15);shg
set(gca,'Ylim',[0 500],'Xlim',[0 5],'Box',true,'FontSize',16);

%% plotting the stripe center (unused for figs). 
% Important bc I botched acquiring the ant and post images
lim = [-0.2 0.25];
xName = 'Fraction of nuclei';

% 15s acqn stripe PDF
figNum = 519;
anal1 = space15anal;
col = colr{1};
ax519 = figure(figNum);
% rep 1
out1 = positionDistribution(anal1,bins);
dist1 = out1.eventBin./(out1.totalBin);dist1(isnan(dist1)) = 0; 
% combined 
handle = bar(bins,mean(dist1,2),'histc');shg
patch('Parent',gca,'Vertices',handle.Vertices,'Faces',handle.Faces,...
    'FaceColor',col,...
    'EdgeColor',colr{1},...
    'CData',handle.CData);
set(gca,'Box',true,'FontSize',16,'XLim',lim,'YLim',[0 1]); 
ylabel(xName)
xlabel('Position relative to stripe center (fraction of AP length)')

% 30s acqn stripe PDF
figNum = 520;
anal1 = spaceAnal1;
anal2 = spaceAnal2;
col = colr{2};
ax520 = figure(figNum);
% rep 1
out1 = positionDistribution(anal1,bins);
dist1 = out1.eventBin./(out1.totalBin);dist1(isnan(dist1)) = 0; 
% rep 2
out2 = positionDistribution(anal2,bins);
dist2 = out2.eventBin./(out2.totalBin);dist2(isnan(dist2)) = 0; 
handle = bar(bins,mean([dist1 dist2],2),'histc');shg
patch('Parent',gca,'Vertices',handle.Vertices,'Faces',handle.Faces,...
    'FaceColor',col,...
    'EdgeColor',colr{2},...
    'CData',handle.CData);
set(gca,'Box',true,'FontSize',16,'XLim',lim,'YLim',[0 1]); 
ylabel(xName)
xlabel('Position relative to stripe center (fraction of AP length)')

% 60s acqn stripe PDF
figNum = 521;
anal1 = space60anal;
col = colr{3};
ax521 = figure(figNum);
% rep 1
out1 = positionDistribution(anal1,bins);
dist1 = out1.eventBin./(out1.totalBin);dist1(isnan(dist1)) = 0; 
% combined 
handle = bar(bins,mean(dist1,2),'histc');shg
patch('Parent',gca,'Vertices',handle.Vertices,'Faces',handle.Faces,...
    'FaceColor',col,...
    'EdgeColor',colr{3},...
    'CData',handle.CData);
set(gca,'Box',true,'FontSize',16,'XLim',lim,'YLim',[0 1]); 
ylabel(xName)
xlabel('Position relative to stripe center (fraction of AP length)')
% 220208 NB: These PDFs look all the same and, importantly, have their
% peaks in the middle two bins, where our "center" nukes live

%% first passage times (SFig 1B):
% plot data first:
fignum = 550;
ax550 = figure(fignum);
repNum = 3;
NtV = zeros(1,repNum);
% 15s
[out,p] = initialFractionBin2(space15anal,bins,boiC,[],15,3000,fignum);
set(p,'Color',colr{1},'LineWidth',1);
hold on
NtV(1) = out.totNukes;
% 30 s
[out,p] = initialFractionBin2(spaceAnal1,spaceAnal2,bins,boiC,[],30,3000,fignum);
set(p,'Color',colr{2},'LineWidth',1);
hold on
NtV(2) = out.totNukes;
% 60 s
[out,p] = initialFractionBin2(space60anal,bins,boiC,[],60,3000,fignum);
set(p,'Color',colr{3},'LineWidth',1);
hold on
NtV(3) = out.totNukes + 1; % an ad hoc hack, as the fit blows up when all the nukes turn on
set(gca,'Box',true,'FontSize',16,'xLim',[0 4000],'yLim',[0 1.01]);

% fit (constraining rate param):
inarg = [4 4 4 0.9 0.9 0.9];
FPrepFit = fminsearch('globalGammaFit6',inarg,[],centerCia,NtV);
% add fits to plot:
xs =0:4000;
rate = 205;
stepsV = FPrepFit(1:3);
AfV = FPrepFit(4:end);
for i = 1:repNum
    gCDF = AfV(i)*gamcdf(xs,stepsV(i),rate);
    hold on;plot(xs,gCDF,'--','LineWidth',0.5,'Color',colr{i},'LineWidth',1);shg
end

%% save
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(ax550,[fp 'fpWTspacePBing22021'],'svg');
%% first bind BS
% generate bootstraps
bs = 10000; 
bsCell = {};
tic;
bsCell{1} = firstEventBootstrapGenerator2(space15anal,bins,boiC,bs);
bsCell{2} = firstEventBootstrapGenerator2(spaceAnal1,spaceAnal2,bins,boiC,bs);
bsCell{3} = firstEventBootstrapGenerator2(space60anal,bins,boiC,bs);
toc
%% BS fit:
inarg = [4 4 4 0.9 0.9 0.9];
tic;
fitErrStruc = firstEventBootstrapFitter4(bsCell,inarg);
toc
% get rid of fit outliers:
fitErrStrucAdjusted = fitErrStruc.fitMat;
fitErrStrucAdjusted(fitErrStrucAdjusted(:,6) > 1,:) = [];
fitErrMeanAdjusted = mean(fitErrStrucAdjusted);
fitErrStdAdjusted = std(fitErrStrucAdjusted);
% 220209 this works and is all set!!

%% saved so you dont ahve to run again:
% save C:\Users\thard\Dropbox\matlab\data\t1p54n.dat nFitErrStruc nBsCell
save C:\Users\tth12\Dropbox\matlab\data\wtSpacePBingFirstPassageBootstrapCell220210.dat bsCell 
%% plot FP bootstraps (SFig 1B)
fignum = 551; 
ax551 = figure(fignum);
% load C:\Users\thard\Dropbox\matlab\data\t1p54n1.dat -mat 
tic;
[~,p] = firstEventBootstrapPlotter2(bsCell{1},space15anal,colr{1},fignum);
hold on;
[~,p] = firstEventBootstrapPlotter2(bsCell{2},spaceAnal,colr{2},fignum);
hold on;
[~,p] = firstEventBootstrapPlotter2(bsCell{3},space60anal,colr{3},fignum);
hold on;
set(gca,'YLim',[0 1.01],'Xlim',[0 4000],'Box',true,'FontSize',16); 
% ylabel('Fraction of nuclei')
% xlabel('Time (s)')
toc;
%% save 
fp = 'C:\Users\tth12\Dropbox\hardenTFfunction2020\figures\figs\';
saveas(ax551,[fp 'fpWTspacePBingBootstraps220210'],'svg');

%other things not used for figures:

%% PDFs of active trxn times:
% initialize
figNum = 560;
averages = zeros(1,3);
% 4 frame/min
% bins = [0,15,30:30:210,270:60:390,500,1000,2000];
bins = [0:30:210,270:60:390,500,1000,2000];
bins = [0:60:240,300:60:420,1000,2000];
bins = [0:60:600,1000,2000];
space15 = pdfPlot(space15cia(:,6),bins,figNum);shg
set(space15.h,'Color',colr{1},'LineWidth',2);
set(space15.j,'Color',colr{1});
averages(1) = mean(space15cia(:,6));
% 2 frame/min
hold on;
bins = [0:30:210,270:60:390,500,1000,2000];
bins = [0:60:240,300:60:420,1000,2000];
bins = [0:60:600,1000,2000];
space30 = pdfPlot(spaceCia(:,6),bins,figNum);shg
set(space30.h,'Color',colr{2},'LineWidth',2);
set(space30.j,'Color',colr{2});
averages(2) = mean(spaceCia(:,6));
% 1 frame/min
bins = [0:60:240,300:60:420,1000,2000];
bins = [0:60:600,1000,2000];
hold on;
space60 = pdfPlot(space60cia(:,6),bins,figNum);shg
set(space60.h,'Color',colr{3},'LineWidth',2);
set(space60.j,'Color',colr{3});
averages(3) = mean(space60cia(:,6));
set(gca,'xLim',[0 1500]);

%% plot median active trxn time v duty cycle:
figure(501);
plot(dutyCycle,averages,'.','MarkerSize',15);shg
set(gca,'Ylim',[0 400],'Xlim',[0 5],'Box',true,'FontSize',16);

%% plot median active trxn time v duty cycle:
cia15dropout  = dropout(space15cia);
averagesDrop(1) = mean(cia15dropout(:,6));
cia30dropout = dropout(spaceCia);
averagesDrop(2) = mean(cia30dropout(:,6));
cia60dropout = dropout(space60cia); 
averagesDrop(3) = mean(cia60dropout(:,6));
figure(502);
plot(dutyCycle,averagesDrop,'.','MarkerSize',15);shg
set(gca,'Ylim',[0 400],'Xlim',[0 5],'Box',true,'FontSize',16);

% ND: looking at frames & not time, it is clear that we are missing short
% lived events with 2 frame/min. I think it besat to make this point with
% spots v t and active trxn times. Get error and look at fit values next
% time. ...think on how to analyze missing short events































