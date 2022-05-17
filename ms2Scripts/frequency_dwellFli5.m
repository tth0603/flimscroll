function [mih,plotHandle] = frequency_dwellFli5(varargin)
% [mih,plotHandle] = frequency_dwellFli2(cia1,cia2,...,timestep,colr,fignum,fitParams[optional],tMx[optional])
%
%This function makes a survival histogram of event frequency by their
%dwell times. At the y-intercept (time zero), the toal frequency
%of events is shown (*basically* read as an on rate). As you go out to
%longer dwell times, events become rarer and the frequency will drop off
%(*basically* an off rate). This kind of plot makes it easy very easy to
%evalaute differences in terms of on and off rates since it is not
%normalized to a total number of events (e.g. if you're trying to show a
%difference between on DNA and off DNA events, the specificity will really
%jump out at you here).
%
% is identical to v3, but omits the time after the last active trxn event
% in summing the idles times (i.e. low_event)
%
%  cia(s) == cumulative interval array of events.
%
%  timestep == determines how often the program will make a step in time to
%check how many events are remaining. If this is smaller then your
%aquisition time, you will get a graph with many plauteaus in it. If it
%is larger, it will have large dropoffs. Best to make this value
%something slightly higher than your actual aqisition time.
%
%  colr == a vector length 3 specifying plot color OR the color code in
%  single quotes eg [0.6 0.6 0.6] OR 'c'
%
%  fignum == the figure number, obvi.
%  tMx == an optional scaler. if we want to put a cap on the length of
%  experiment time to compare between different experiments where the
%  observation times during NC14 differ. 
%
%version 3.2 (Nov 2015)
%updated Harden 2018 for flimscroll data
                                                            
%Defining vars
numOfReps = length(varargin) - 5;
timestep = varargin{numOfReps + 1};
colr = varargin{numOfReps + 2};
fignum = varargin{numOfReps + 3};
sints=[]; 
low_events = [];
events = [];
expLength = zeros(numOfReps,1);
%find the stuff for each replicate:
for i = 1:numOfReps
    % get cia, rid of dropouts:
    cia = varargin{i};
    % get the experiment length of each cia:
    lastEvents = cia(cia(:,2) == 3,:);
    lastTimeV = lastEvents(:,7) + lastEvents(:,6); % time they begin rel to NC 14 start plus the event length
    lastTime = max(lastTimeV);
%     expLength(i) = lastTime; % for TSing
    % get the time of the first event as the start of the experiment
    firstTimeV = cia(cia(:,2) == -3,7); %a vector
    % pick out the high event durations (easy):
    events = [events; cia(:,6)];
    % get the low event durations:
    nukes = unique(cia(:,1));
    len = length(nukes);
    for j = 1:len
        mat2 = cia(cia(:,1) == nukes(j),:);
        % get the inter events:
%         intervals = diff(mat2(:,7)) - mat2(1:end-1,6) - 60; % start time of next event minus start time of this event plus this even't duration. The -60 is to get rid of the lag at the start of the survival curve
        intervals = diff(mat2(:,7)) - mat2(1:end-1,6); % get rid of the -60 for now
        % get the low event at the start:
    %     strtInt = mat2(1,7) - firstTime; %time of this nukes first event minus the time we first detect a spot in this experiment
        % get the low event at the end:
        if ~isempty(varargin{numOfReps + 5}) % use a universal last time:
            tMx = varargin{numOfReps + 5};
            if mat2(end,7) + mat2(end,6) < tMx
                endInt = tMx - (mat2(end,7) + mat2(end,6));
            else 
                endInt = [];
            end
        else % use the end of the last event:
            if lastTime < 4000
                endInt = lastTime - (mat2(end,7) + mat2(end,6)); % when this experiment's last event ends minus (this nukes last evnt strt time + the length of that event)
            else % this is an adjustment to correct for one long WT experiment. for this one case, we change the length of the expt to be the average of all 10 experiments used. be careful here
                lastTime = 3.2730e+03;
                endInt = lastTime - (mat2(end,7) + mat2(end,6));
            end
        end
        % or, omit these events:
%         endInt = [];
%         low_events = [low_events; intervals; endInt]; % use this if you wanna include the end events as in v3
        low_events = [low_events; intervals];
    end
end
% and now combing both cias:
low_times = sum(low_events);                                               %...and suming their times

for i=0:timestep:max(events)                                            %A for loop for making a survival histogram. Timesteps are input by user.
logik=events>=i;                                                           %Counts events whose dwell time is greater than or equal to the current timestep...
sints=[sints;i (sum(logik))];                                              %...and makes an array of those dwell times and those events called 'sints'
end

sints = [sints(:,1), sints(:,2)/low_times];                                %Changes total number of events to frequencies

%Below plots the survival data on a linear axis
figure(fignum);
hold on;
% plot points:
% plotHandle = plot(sints(:,1),(sints(:,2)),'MarkerFaceColor',colr,'MarkerSize',3,'Marker','o',...
%     'LineStyle','none',...
%     'Color',colr);
% plot stairs:
plotHandle = stairs(sints(:,1),sints(:,2),'LineWidth',1,'Color',colr);

set(gca,'Box',true,'FontSize',16);
set(gca,'YScale','log')
% xlabel('Dwell times(s)');ylabel('Frequency(s ^-^1)');
%Below are common ouputs that I find useful to have
mih.sints = sints;
totalFrequency1 = length(events)/low_times; % this counts the first event in the freq calc, which is fine for the plot, but not in reporting idel periods, as we already account for the first event in the first passage plot
totalFrequency2 = length(low_events)/low_times; % This looks at frequency following the end of the first event. Which is appropriate for reporting idle times. 
mih.totalFrequency = totalFrequency2;
mih.eventNumber = length(events);
mih.lowTimes = low_times;
mih.expLength = expLength;

if ~isempty(varargin{numOfReps + 4})
    %from survivalPLotFromVector.m with som revisions
    inarg = varargin{numOfReps + 4};
    tm = 1;tx = 3000;
    params=fminsearch('expfalltwo_mxl',inarg,[],events,tm,tx); 
    corr=expfalltwo_mxl_correction(params,events,tm,tx);      %correction factors for the fact that want to characterize the entire phenomena, but only observe it for a part of the time
        %'corr' is a struc: 
        % corr.CorrectionFactor       % = (  1/( a*A + (1-a)*B ) 
        % corr.Nzero                  % = actual # of events in data set (observed+not observed)  
        % corr.Offset            % # remaining events at time = tx 
    ap=params(1);
    a=1/(1+ap^2);
    k1=params(2);
    k2=params(3);  %defines rate constants for fit
    tfit=[0:max(events)];
    SPfit=corr.CorrectionFactor*(a*exp(-tfit*k1)+(1-a)*exp(-tfit*k2));  %This is the survival 
%     SPfit=(a*exp(-tfit*k1)+(1-a)*exp(-tfit*k2));
%this is the expn that we actually use to fit to account for exp
%resolution:
%     A=( exp(-tm*k1) - exp(-tx*k1) );
%     B = ( exp(-tm*k2) - exp(-tx*k2) );
%     SPfit=(  1/( a*A + (1-a)*B  ) )*...
%                                  ( a*k1*exp(-tfit*k1)+(1-a)*k2*exp(-tfit*k2) )+...
%                                  0;
%     SPfit = a*k1*exp(-tfit*k1)+(1-a)*k2*exp(-tfit*k2);
    hold on;
%     plotHandle = plot(tfit,SPfit*totalFrequency,'Color',[0.4 0.4 0.4]);shg
    plotHandle = plot(tfit,SPfit*totalFrequency1,'--','Color',colr,'LineWidth',1);shg
    mih.k1=params(2);
    mih.k2=params(3);
    mih.a = a;
end

end