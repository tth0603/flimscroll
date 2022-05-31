function [output] = frequency_dwell(cia,timestep,color_shape,fignum)
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
%  color_symbol == the one letter color code followed by the symbol. Needs 
%to be entered with '' around the argument, eg 'r-x' or 'ko'.
%
%  fignum == the figure number, obvi.
%
%version 3.2 (Nov 2015)
sints=[];                                                                  %Defining vars

events = cia((cia(:,1)== -3 | cia(:,1)== 1 | cia(:,1)== 3),5);             %Picking out the high events    
low_events = cia((cia(:,1) == -2 | cia(:,1) == 0 | cia(:,1) == 2),5);      %Picking out low events...
low_times = sum(low_events);                                               %...and suming their times

for i=0:timestep:max(events)                                               %A for loop for making a survival histogram. Timesteps are input by user.
logik=events>=i;                                                           %Counts events whose dwell time is greater than or equal to the current timestep...
sints=[sints;i (sum(logik))];                                              %...and makes an array of those dwell times and those events called 'sints'
end

sints = [sints(:,1), sints(:,2)/low_times];                                %Changes total number of events to frequencies

%Below plots the survival data on a linear axis
hold on;figure(fignum);hold on;plot(sints(:,1),(sints(:,2)),color_shape);xlabel('Dwell times(s)');ylabel('Frequency(s ^-^1)');title('Survival histogram of dwell times')

%Below are common ouputs that I find useful to have
output.sints = sints;
output.total_frequency = length(events)/low_times;

end
