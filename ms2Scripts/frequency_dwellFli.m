function [mih,plotHandle] = frequency_dwellFli(cia,timestep,maxTime,colr,fignum)
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
%  cia == cumulative interval array of events.
%
%  timestep == determines how often the program will make a step in time to
%check how many events are remaining. If this is smaller then your
%aquisition time, you will get a graph with many plauteaus in it. If it
%is larger, it will have large dropoffs. Best to make this value
%something slightly higher than your actual aqisition time.
%
%  maxTime == the length of the experiment in seconds (a scaler). if empty, then
%  default to when the last event in cia ends. if comparing multiple
%  experiemnts, use the time of the longest experiment.
%
%  colr == a vector length 3 specifying plot color OR the color code in
%  single quotes eg [0.6 0.6 0.6] OR 'c'
%
%  fignum == the figure number, obvi.
%
%version 3.2 (Nov 2015)
%updated Harden 2018 for flimscroll data
sints=[];                                                                  %Defining vars
% get rid of dropouts:
cia = dropout(cia);

%TH: if I get this straight, we want to pick ou the on times, then sum the
%total off times, starting not from the start of nc 14 but instead at teh 
%appearance of teh first spot to account for the non-single exp nature of
%the first on event. The latter requires findind the time of the last acqn
%frame. I guess we'll get this from a nukeMat. so the beginning of each
%'recording' will be different for each nuke, but the end of the recording
%will be the same? does this make sense what with the weirdness at the end
%of nc14? Right now I am not standardizing when I 'turn off' the recording
%from experiment to experiment. The best would be to use the time when we
%see the last spot fo away. but damn I think it'll be hard to interpret
%these results.
% do that here before you forget:

% get the time of the first event as the start of the experiment
firstTimes = cia(cia(:,2) == -3,7); %a vector
%make the first event occur at t = 0;
firstTime = min(firstTimes); %harden 202009: you can't have your first event start at zero 
% get the time that the last event ends (a bit trickier)\:
if isempty(maxTime)
    lastEvents = cia(cia(:,2) == 3,:);
    lastTimes = lastEvents(:,7) + lastEvents(:,6); % time they begin rel to NC 14 start plus the event length
    lastTime = max(lastTimes);
else 
    lastTime = maxTime;
end
% lastTime = 3.6150e+03;
% pick out the high event durations (easy):
events = cia(:,6);
% get the low event durations:
nukes = unique(cia(:,1));
len = length(nukes);
low_events = [];
for j = 1:len
    mat2 = cia(cia(:,1) == nukes(j),:);
    % get the inter events:
    intervals = diff(mat2(:,7)) - mat2(1:end-1,6) - 60; % start time of next event minus start time of this event plus this even't duration. The -60 is to get rid of the lag at the start of the survival curve
%     intervals = diff(mat2(:,7)) - mat2(1:end-1,6); % get rid of the -60 for now
    % get the low event at the start:
%     strtInt = mat2(1,7) - firstTime; %time of this nukes first event minus the time we first detect a spot in this experiment
    strtInt = [];
    % get the low event at the end:
    endInt = lastTime - (mat2(end,7) + mat2(end,6)); % when this experiment's last event ends minus (this nukes last evnt strt time + the length of that event)
%     endInt = [];
    low_events = [low_events; intervals; strtInt; endInt];
end
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
xlabel('Dwell times(s)');ylabel('Frequency(s ^-^1)');
%Below are common ouputs that I find useful to have
mih.sints = sints;
mih.totalFrequency = length(events)/low_times;

end