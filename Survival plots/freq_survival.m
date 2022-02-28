function [freq_survival] =freq_survival(cia,control_cia)
a = cia(cia(:,1)==1,[1:7]);
j = control_cia(control_cia(:,1)==1,[1:7]);
k=j(:,5);
d=a(:,5);
if (max(k)>=max(d))
    max_time=max(k);
else
    max_time=max(d);
end
e=length(d);
f = cia(cia(:,1) == -2 | cia(:,1) == 0 | cia(:,1) == 2,[5]);
g=sum(f); %total low time
h=e/g; %frequency of events 

sints=[]; %below will make a survival plot
for i=0:.1:max_time
logik=d(:,1)>i;
sints=[sints;i (sum(logik)/g)];
end

j = control_cia(control_cia(:,1)==1,[1:7]);
k=j(:,5);
l=length(k);
m = control_cia(control_cia(:,1) == -2 | control_cia(:,1) == 0 | control_cia(:,1) == 2,[5]);
n=sum(m); %total low time
o=l/n; %frequency of events 

sints_control=[]; %below will make a survival plot
for i=0:.1:max_time
logik=k(:,1)>i;
sints_control=[sints_control;i (sum(logik)/n)];
end
y=length(a(:,1));
z=length(j(:,1));

freq_survival = struct('frequency_of_cia_events',{(1/h)},'number_of_cia_events',{y},'number_of_cia_AOIs',{max(cia(:,7))},'survival_cia',{sints},'frequency_of_control_events',{1/o},'number_of_control_events',{z},'number_of_control_AOIs',{max(control_cia(:,7))},'survival_control',{sints_control})
figure(4);plot(sints_control(:,1),(sints_control(:,2)),'-b');hold on; plot(sints(:,1),(sints(:,2)),'-g');xlabel('Dwell Time (s)');ylabel('Frequency of Remaining Events (1/s)');title('Survival Plot of Event Frequency')