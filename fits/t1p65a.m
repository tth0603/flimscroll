%% fitting these distributions to all data
% so this fitting actually works, but is highly dependent on initial
% cnditions. As such we'll use the fiddly fits to determine the initial
% guesses
% paramMat = [0.28 0.0025 0.0010;...
%     1.00 0.015 0.0023;...
%     1.2 0.0035 0.0031;...
%     0.475 0.0018 0.0012;...
%     0.40 0.002 0.0013];
paramMat = [0.28 0.0025 0.001;... %it'll have to do. There really isnat any fitting going on [0.2812    0.0026    0.0010]
    0.93 0.015 0.0023;... %good; result: [0.9893    0.0148    0.0023]
    0.53 0.0035 0.0031;... % PG result: [1.22 0.0023 0.0052]
    0.5 0.002 0.001;... %PG, but not a lot of fitting going on result: [0.5027 0.002 0.001]
    0.40 0.002 0.0013]; %good. result [0.9544    0.0011    0.0019]
options = optimset('MaxIter',10^7,'MaxFunEvals',10^7); %for fits
fitMat = zeros(5,3);
tic;
for i = 1:5
    dataSet = i; %1 - 5.. [1.neutral 2.wt 3.zld 4.bcd 5.dst]
    cia = ciaCell{dataSet};
    first = cia(cia(:,2) == -3,7); %a vector
    %make the first event occur at t = 0;
    first = first - min(first);
    params = paramMat(i,:);
%     [fitParms,fval,ef,out] = fminsearch('twoStepBindFitV1',params,options,first);
    [fitParms] = fminsearch('twoStepBindFitV1',params,options,first);
    fitMat(i,:) = fitParms;
end
toc %takes about 13 minutes on my laptop to fit all 5, ~ 37s of which is for wt

%result:
%     0.2812    0.0026    0.0010
%     0.9893    0.0148    0.0023
%     1.2220    0.0023    0.0052
%     0.5027    0.0020    0.0010
%     0.9544    0.0011    0.0019

%% plot fits with curves
xs = 1:4000;
fignum = 5;
% for i = dataSet
for i = 1:5
    [~,p] = initialFraction2(ciaCell{i},238,[],fignum + 1); % 238 is # of Wt events
    set(p,'Color',colr{i},'LineWidth',1);
    fA = fitMat(i,1);fk0 = fitMat(i,2);fk1 = fitMat(i,3);
    fit = fA*( fk0*fk1/(fk0 - fk1).*( (1/fk0)*exp(-fk0.*xs) - (1/fk1)*exp(-fk1.*xs) ) + 1 );
    hold on;p2 = plot(xs,fit,'.');
    set(p2,'Color',colr{i});   
end
% set(gca,'YLim',[0 0.9],'Box',true,'FontSize',16); % pair with "A"
ylabel('fraction of total nuclei in the "pattern"')
xlabel('time (s)')

%% figuring out individual fits
indx =  1; %pick which data set to use, neutral, wt, zld, bcd or dst
paramMat = [0.28 0.0025 0.001;... %it'll have to do. There really isnat any fitting going on [0.2812    0.0026    0.0010]
    0.93 0.015 0.0023;... %good; result: [0.9893    0.0148    0.0023]
    0.53 0.0035 0.0031;... % PG result: [1.22 0.0023 0.0052]
    0.5 0.002 0.001;... %PG, but not a lot of fitting going on result: [0.5027 0.002 0.001]
    0.40 0.002 0.0013]; %good. result [0.9544    0.0011    0.0019]
options = optimset('MaxIter',10^7,'MaxFunEvals',10^7); %for fits
fitMat = zeros(5,3);
tic;
for i = indx
    dataSet = i; %1 - 5.. [1.neutral 2.wt 3.zld 4.bcd 5.dst]
    cia = ciaCell{dataSet};
    first = cia(cia(:,2) == -3,7); %a vector
    %make the first event occur at t = 0;
    first = first - min(first);
    params = paramMat(i,:);
    [fitParms,fval,ef,out] = fminsearch('twoStepBindFitV1',params,options,first);
    fitMat(i,:) = fitParms;
end
toc

%% plot result
fignum = 5;
% for i = dataSet
for i = indx
    [~,p] = initialFraction2(ciaCell{i},238,[],fignum + 1); % 238 is # of Wt events
    set(p,'Color',colr{i},'LineWidth',1);
    fA = fitMat(i,1);fk0 = fitMat(i,2);fk1 = fitMat(i,3);
    fit = fA*( fk0*fk1/(fk0 - fk1).*( (1/fk0)*exp(-fk0.*xs) - (1/fk1)*exp(-fk1.*xs) ) + 1 );
    hold on;p2 = plot(xs,fit,'.');
    set(p2,'Color',colr{i});   
end
% set(gca,'YLim',[0 0.9],'Box',true,'FontSize',16); % pair with "A"
ylabel('fraction of total nuclei in the "pattern"')
xlabel('time (s)')

%% fiddle fit:
% this is to play around 
k0 = 0.0023;
k1 = 0.0052;
A = 1.22;
ys = A*(k0*k1/(k0 - k1).*(exp(-k1.*xs) - exp(-k0.*xs))); % just tack on a coefficient. can fiddle with formalities later
hold on;
plot(xs,cumsum(ys));shg %we integrate via cumsum
%keepers
% neutral: 0.28 0.0025 0.0010 or 0.30 0.001 0.002 
% wt: 1.00 0.015 0.0023 
% Zld: 1.2 0.0035 0.0031 
% Bcd: 0.475 0.0018 0.0012 or 0.475 0.001 0.002 
% Dst: 0.87 0.0020 0.0013 

















