    %playing with gamma functions
% you'll have to load the workspace from t1p54b or t1p54e to run this
%% lets just try and plot some:
tst = gammainc(2,xs/tau);
figure(1);plot(xs,1./tst);shg
%% from lagha 2018
xs = 0:4000;
tau = 500; %a guess at teh time const
% probabilities:
p1 = 1; p2 = 1; p3 = 1;
%% one and two step:
tau = 1000;
p1 = 0.9; p2 = 0.9;
g1 = (1 - exp(-xs./tau));
g2 = p1.*(1 - exp(-xs./tau)) + p2.*(1/gamma(2)).*gammainc(xs/tau,2); %not dure why this isn't a cdf, as in it doesnt asymtote to 1
% figure(2);plot(xs,g2);shg
hold on;plot(xs,g2);shg
% hold on; plot(xs,g1);
%% 3 step:
tau = 100;
p1 = 0.33; p2 = 0.33; p3 = 0.33;
g3 = p1.*(1 - exp(-xs./tau)) + p2.*(1/gamma(2)).*gammainc(xs/tau,2) + p3.*(1/gamma(3)).*gammainc(xs/tau,3);
% figure(3);plot(xs,g3);shg
hold on;plot(xs,g3);shg
% figure(3);plot(xs,cumsum(g3));shg

%% play with incomplete gamma fn
tau = 500;
g3inc = gammainc(xs/tau,3);
figure(4);plot(xs,g3inc);shg

%% use ML fitting to find A number of steps and B time constant:
%make the events start at origin, make mat for fit values
gFitMat = zeros(5,2);
tic;
for i = 1:5
    dataSet = i; %1 - 5.. [1.neutral 2.wt 3.zld 4.bcd 5.dst]
    cia = ciaCell{dataSet};
    first = cia(cia(:,2) == -3,7); %a vector
    %make the first event occur at t = 0;
%     first = first - min(first);
    % if you want to get rid of zeros:
    first(first == 0) = 1;
    [gFitMat(i,:) ci{i}] = gamfit(first,0.1); %result of each iteration is: [a b] gamma parameters for shape parameter a (number of steps) and scale parameter b (rate constant)
end
toc

%% now compare the fits with the curves:
fignum = 12;
%get the normalization:
for i = 1:5
    dataCells{i}.Nt = size(ciaCell{i}(ciaCell{i}(:,2) == -3,:),1); 
end
% plot data:
for i = 1:5
    [~,p] = initialFraction2(ciaCell{i},dataCells{i}.Nt,[],fignum); % "C". all norm'd to one, all start at origin
    [~,p] = initialFraction(ciaCell{i},dataCells{i}.Nt,[],fignum); % "C". all norm'd to one, all start times
    set(p,'Color',colr{i},'LineWidth',1);
    hold on
end
set(gca,'YLim',[0 1.01],'Xlim',[0 4000],'Box',true,'FontSize',16); 
ylabel('fraction of total nuclei in the "pattern"')
xlabel('time (s)')
% plot gamma fit
xs =0:4000;
for i = 1:5
    gCDF = gamcdf(xs,gFitMat(i,1),gFitMat(i,2));
    hold on;plot(xs,gCDF,'Color',colr{i},'LineWidth',1);shg
end
% NB: this fitting works pretty well, but it normalizes everything to one.
% try normalizing them to fraction of nukes within the nukes.

% %the following ins from t1p54k. you need to load that workspace to run is
% 
% none of this commented out stuff made it into the paper, but is playing
% around with first apssage times:
% %% initial ts - two equal steps
% figNum = 36;
% fitInputs = [0.1 100]; 
% % neutral
% anal1 = nAnal1;
% anal2 = nAnal2;
% [ini,wtH] = initialFractionBin(anal1,anal2,bins,boiC,fitInputs,30,3000,figNum);
% set(wtH,'Color',colr{1});
% ini
% % space
% anal1 = spaceAnal1;
% anal2 = spaceAnal2;
% [ini,wtH] = initialFractionBin(anal1,anal2,bins,boiC,fitInputs,30,3000,figNum);
% set(wtH,'Color',colr{2});
% ini
% % zld
% anal1 = zldAnal1;
% anal2 = zldAnal2;
% [ini,wtH] = initialFractionBin(anal1,anal2,bins,boiC,fitInputs,30,3000,figNum);
% set(wtH,'Color',colr{3});
% ini
% % bcd
% anal1 = bcdAnal1;
% anal2 = bcdAnal2;
% [ini,wtH] = initialFractionBin(anal1,anal2,bins,boiC,fitInputs,30,3000,figNum);
% set(wtH,'Color',colr{4});
% ini
% % dst
% anal1 = dAnal1;
% anal2 = dAnal2;
% [ini,wtH] = initialFractionBin(anal1,anal2,bins,boiC,fitInputs,30,3000,figNum);
% set(wtH,'Color',colr{5});
% ini
% 
% set(gca,'Box',true,'FontSize',16,'XLim',[0 3500],'YLim',[0 1]);
% ylabel('Fraction of nuclei in stripe domain')
% xlabel('time (s)')

%% gamma fns - all nukes
% %make the events start at origin, make mat for fit values
% gFitMat = zeros(5,2);
% tic;
% for i = 1:5
%     cia = ciaCell{i};
%     first = cia(cia(:,2) == -3,7); %a vector
%     %make the first event occur at t = 0;
%     first = first - min(first);
%     % if you want to get rid of zeros:
%     first(first == 0) = 1;
%     [gFitMat(i,:) ci{i}] = gamfit(first,0.1); %result of each iteration is: [a b] gamma parameters for shape parameter a (number of steps) and scale parameter b (rate constant)
% end
% toc
% % now compare the fits with the curves:
% fignum = 12;
% %get the normalization:
% for i = 1:5
%     dataCells{i}.Nt = size(ciaCell{i}(ciaCell{i}(:,2) == -3,:),1); 
% end
% % plot data:
% for i = 1:5
%     [~,p] = initialFraction2(ciaCell{i},dataCells{i}.Nt,[],fignum); % "C". all norm'd to one, all start at origin
%     set(p,'Color',colr{i},'LineWidth',1);
%     hold on
% end
% set(gca,'YLim',[0 1.01],'Xlim',[0 4000],'Box',true,'FontSize',16); 
% ylabel('fraction of total nuclei in the "pattern"')
% xlabel('time (s)')
% % plot gamma fit
% xs =0:4000;
% for i = 1:5
%     gCDF = gamcdf(xs,gFitMat(i,1),gFitMat(i,2));
%     hold on;plot(xs,gCDF,'Color',colr{i},'LineWidth',1);shg
% end
% 
% %% gamma fns - all nukes - all times
% %make the events start at origin, make mat for fit values
% fignum = 13;
% gFitMat = zeros(5,2);
% tic;
% for i = 1:5
%     cia = ciaCell{i};
%     first = cia(cia(:,2) == -3,7); %a vector
%     [gFitMat(i,:) ci{i}] = gamfit(first,0.1); %result of each iteration is: [a b] gamma parameters for shape parameter a (number of steps) and scale parameter b (rate constant)
% end
% toc
% % now compare the fits with the curves:
% %get the normalization:
% for i = 1:5
%     dataCells{i}.Nt = size(ciaCell{i}(ciaCell{i}(:,2) == -3,:),1); 
% end
% % plot data:
% for i = 1:5
%     [~,p] = initialFraction(ciaCell{i},dataCells{i}.Nt,[],fignum); % "C". all norm'd to one, all start at origin
%     set(p,'Color',colr{i},'LineWidth',1);
%     hold on
% end
% set(gca,'YLim',[0 1.01],'Xlim',[0 4000],'Box',true,'FontSize',16); 
% ylabel('fraction of total nuclei in the "pattern"')
% xlabel('time (s)')
% % plot gamma fit
% xs =0:4000;
% for i = 1:5
%     gCDF = gamcdf(xs,gFitMat(i,1),gFitMat(i,2));
%     hold on;plot(xs,gCDF,'Color',colr{i},'LineWidth',1);shg
% end
% 
% %% inital ts - gamma fns - center nukes - all times - norm to one
% %make the events start at origin, norm'd to one.
% fignum = 14;
% gFitMat = zeros(5,2);
% tic;
% for i = 1:5
%     cia = centerCia{i};
%     first = cia(cia(:,2) == -3,7); %a vector
%     [gFitMat(i,:) ci{i}] = gamfit(first,0.1); %result of each iteration is: [a b] gamma parameters for shape parameter a (number of steps) and scale parameter b (rate constant)
% end
% toc
% 
% %get the normalization:
% for i = 1:5
%     dataCells{i}.Nt = size(centerCia{i}(centerCia{i}(:,2) == -3,:),1); 
% end
% % plot data:
% for i = 1:5
%     [~,p] = initialFraction(centerCia{i},dataCells{i}.Nt,[],fignum); % "C". all norm'd to one, all start at origin
%     set(p,'Color',colr{i},'LineWidth',1);
%     hold on
% end
% set(gca,'YLim',[0 1.01],'Xlim',[0 4000],'Box',true,'FontSize',16); 
% ylabel('fraction of total nuclei in the "pattern"')
% xlabel('time (s)')
% % plot gamma fit
% xs =0:4000;
% for i = 1:5
%     gCDF = gamcdf(xs,gFitMat(i,1),gFitMat(i,2));
%     hold on;plot(xs,gCDF,'Color',k,'LineWidth',0.5);shg
% end
% 
% %% inital ts - gamma fns - center nukes - norm'd to nukes/bin
% %make the events start at origin, make mat for fit values
% fignum = 14;
% gFitMat = zeros(5,2);
% tic;
% for i = 1:5
%     cia = centerCia{i};
%     first = cia(cia(:,2) == -3,7); %a vector
% %     %make the first event occur at t = 0;
% %     first = first - min(first);
% %     % if you want to get rid of zeros:
% %     first(first == 0) = 1;
%     [gFitMat(i,:) ci{i}] = gamfit(first,0.1); %result of each iteration is: [a b] gamma parameters for shape parameter a (number of steps) and scale parameter b (rate constant)
% end
% toc
% 
% % plot data:
% for i = 1:5
%     [~,p] = initialFractionBin2(analCell1{i},analCell2{i},bins,boiC,[],30,3000,fignum);
%     set(p,'Color',colr{i},'LineWidth',1);
%     hold on
% end
% set(gca,'YLim',[0 1.01],'Xlim',[0 4000],'Box',true,'FontSize',16); 
% ylabel('fraction of total nuclei in the "pattern"')
% xlabel('time (s)')
% Af = [0.4395    0.8977    0.8600    0.7442    0.9295]; % these are from fitting the initial curves with initialFractionBin:
% % fitInputs = [0.1 100]; 
% % Af = zeros(1,5);
% % for i = 1:5
% %     [ini,p] = initialFractionBin(analCell1{i},analCell2{i},bins,boiC,fitInputs,30,3000,fignum+1);
% %     Af(i) = ini.Af;
% %     set(p,'Color',colr{i},'LineWidth',1);
% %     hold on
% % end
% % plot gamma fit
% xs =0:4000;
% for i = 1:5
%     gCDF = gamcdf(xs,gFitMat(i,1),gFitMat(i,2));
%     hold on;plot(xs,Af(i).*gCDF,'Color',colr{i},'LineWidth',1);shg
% end
% 
% %% global gamma fit: global rate - norm'd to nukes/bin
% % center nukes, norm'd to one, Fit the rate (scale
% % parameter,B) constant globally:
% inarg = [100 1 1 1 1 1]; % [globalRate step1 step2 ...];
% gFit = fminsearch('globalGammaFit1',inarg,[],centerCia)
% %plot the curves:
% fignum = 15;
% for i = 1:5
%     dataCells{i}.Nt = size(centerCia{i}(centerCia{i}(:,2) == -3,:),1); 
% end
% % plot data:
% for i = 1:5
%     [~,p] = initialFractionBin2(analCell1{i},analCell2{i},bins,boiC,[],30,3000,fignum);
%     set(p,'Color',colr{i},'LineWidth',1);
%     hold on
% end
% set(gca,'YLim',[0 1.01],'Xlim',[0 4000],'Box',true,'FontSize',16); 
% ylabel('Fraction of nuclei')
% xlabel('Time (s)')
% % and fits:
% Af = [0.4395    0.8977    0.8600    0.7442    0.9295];
% xs =0:4000;
% for i = 1:5
%     gCDF = gamcdf(xs,gFit(i + 1),gFit(1));
%     hold on;plot(xs,Af(i).*gCDF,'Color',colr{i},'LineWidth',1);shg
% end
% 
% %% global gamma fit: global step
% % center nukes, norm'd to one, Fit the no of steps (shape
% % parameter,A) globally:
% inarg = [1 100 100 100 100 100]; % [globalSte rate1 rate2...];
% gFit2 = fminsearch('globalGammaFit2',inarg,[],centerCia);
% %plot the curves:
% fignum = 16;
% for i = 1:5
%     dataCells{i}.Nt = size(centerCia{i}(centerCia{i}(:,2) == -3,:),1); 
% end
% % plot data:
% for i = 1:5
%     [~,p] = initialFraction(centerCia{i},dataCells{i}.Nt,[],fignum); % "C". all norm'd to one, all start at origin
%     set(p,'Color',colr{i},'LineWidth',1);
%     hold on
% end
% set(gca,'YLim',[0 1.01],'Xlim',[0 4000],'Box',true,'FontSize',16); 
% ylabel('fraction of total nuclei in the "pattern"')
% xlabel('time (s)')
% % and fits:
% xs =0:4000;
% for i = 1:5
%     gCDF = gamcdf(xs,gFit2(1),gFit2(i + 1));
%     hold on;plot(xs,gCDF,'Color',colr{i},'LineWidth',1);shg
% end
% % NB: holding the rate const fits better. I am not sure which
% % makes more sense biologically. perhaps holding step const does, but that
% % fit is bad enough (esp for neutral) that I don't think we can choose that
% % one. 
% 
% %% global gamma fit: global rate, global postion
% % center nukes, norm'd to one, Fit the rate (scale
% % parameter,B) and location mu globally:
% inarg = [100 1 1 1 1 1 10]; % [globalRate step1 step2 ... mu];
% gFit3 = fminsearch('globalGammaFit3',inarg,[],centerCia)
% % mn = 0.01; mx = 20;
% % gFit3 = fminsearchbnd('globalGammaFit3',inarg,[1 mn mn mn mn mn 0],[3000 mx mx mx mx mx 1000],[],centerCia);
% %plot the curves:
% fignum = 17;
% for i = 1:5
%     dataCells{i}.Nt = size(centerCia{i}(centerCia{i}(:,2) == -3,:),1); 
% end
% % plot data:
% for i = 1:5
%     [~,p] = initialFraction(centerCia{i},dataCells{i}.Nt,[],fignum); % "C". all norm'd to one, all start at origin
%     set(p,'Color',colr{i},'LineWidth',1);
%     hold on
% end
% % set(gca,'YLim',[0 1.01],'Xlim',[0 4000],'Box',true,'FontSize',16); 
% ylabel('fraction of total nuclei in the "pattern"')
% xlabel('time (s)')
% % and fits:
% xs =0:4000;
% for i = 1:5
%     % A (gamma) is shape (no of steps), B (beta) is scale (rate), mu is location
%     mu = abs(gFit3(end));
%     B = abs(gFit3(1));
%     A = abs(gFit3(2:end - 1));
%     gCDF = cumsum( ( ((xs - mu)./B).^(A(i) - 1).*exp( -(xs - mu)./B ) )./(B*gamma( A(i) )) );
%     hold on;plot(xs,gCDF,'Color',colr{i},'LineWidth',1);shg
% end
% % NB: while this method doesn't work perfectly, it works well enough that I
% % am confident in saying that if there is any location shift in the data it
% % is less than 5s. 













































