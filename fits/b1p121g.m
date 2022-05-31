load C:/matlab/data/b1p121c.dat -mat %Interval Data File, oligo spot duration, 1 EXP fit


%The concentration for this experiment is [RNAP] = 0.8 nM, the 
%rate constant that tis file spits out is a second order rate constant
%and thereby must be scaled by the RNAP conc.

cia=Intervals.CumulativeIntervalArray; %inspect Intervals variable using 
                                       %Intervals.CumulativeIntervalArrayDescription
logik=(cia(:,1)==0)|(cia(:,1)==1);

data=[cia(logik,5) cia(logik,7) cia(logik,1)];

L=length(data);

ints=[];
for i=1:L-1;
    if (data(i,2)==data(i+1,2)) && (data(i,3)==1);
    ints(i)=data(i,1)+data(i+1,1);
    end
end
ints=ints(ints~=0);

intervals=ints';

%disp('Sample Intervals:');
%disp(intervals(1:10));  %Defines the intervals data to be plotted

%lg=intervals<200;
%intervals=intervals(lg);         
                            %removes two outliers above 300.  May want to 
                            %play with this. Changes the lenth of
                            %'sints' from 290 to 288.

ap=1;
a=1/(1+ap^2);
tau=0.001;    %tau refers to rate constants, change to time constants by    
             %edit expfalltwo_mxl       
tm=0.5;
tx=max(intervals);      %tm is (1/2)*frame length, tm is the longest event 
                      %in the data set

params=fminsearch('expfallone_mxl',tau,[],intervals,tm,tx);

k=params(1);  %gives rate constants for fit
disp('Rate Constant (1/s):');
disp(k);
t=1./k;
disp('Time Constant (s):');
disp(t);

%sints=sort(intervals);  %this sorts the intervals data by length

tfit=0:0.5:max(intervals);

SPfit=exp(-tfit*k);  %This is the survival curve model
                                                     
SPx=exp(-tx*k);  
SPm=exp(-tm*k);  %survival probability evaluated at tm and tx

%for i=1:length(sints)
    %SPmut(i)=(1-i/length(sints))*(SPm-SPx)+SPx;    
%end
                                       %This is what the data set will be 
                                      %plotted against.  Instead of 
                                      %straight time, it is plotted against 
                                      %the time intervals of events, then 
                                      %scaled to account for the fact that
                                      %we are describing the entire
                                      %phenomonon, but only observing a
                                     %small window of time 
                                     
sints=[];
for i=1:max(intervals)
logik=intervals>i;
sints=[sints;i sum(logik)];
end

l=max(sints(:,2));

figure(70);plot(tfit,log(SPfit),sints(:,1),log(sints(:,2)/l),'r');shg

%figure(6);plot(tfit,log(SPfit),sints,log(SPmut));shg