%This will give initial binding curves (exp decay) for both the WT and the mutant landings.
%must have interval data file for each.  Start by loading the interval file:

load t:\matlab\data\b1p63a.dat -mat   %Creates 'Intervals' variable.

cia=Intervals.CumulativeIntervalArray;
logik=cia(:,1)==-2;
ints=cia(logik,4);                  %These are the times before initial binding

%disp('Sample WT intervals:');
%disp(ints(1:10));

ibind=[];                        %creates the variable that will be the 
                                %sorted initial binding times
for indx=0:1690
    logik=ints>indx;
    ibind=[ibind;indx sum(logik)];
end                             %Creats a vector of the sorted times
                                %Differs from the sorted intervals that Ben
                                %showed you in b1p64a.m.  This vector gives
                                %at each time interval (every second) gives
                                %the number of landings that are of that
                                %duration or longer
%Repeat for Mutant Data:
load t:\matlab\data\b1p62c.dat -mat

ciamut=Intervals.CumulativeIntervalArray;
logik=ciamut(:,1)==-2;
intsmut=ciamut(logik,4);

%disp('Sample MUT intervals:');
%disp(intsmut(1:10));

ibindmut=[];

for indx=0:1690
    logik=intsmut>indx;
    ibindmut=[ibindmut;indx sum(logik)];
end  

figure(5);plot(ibind(:,1),log(ibind(:,2)/max(ibind(:,2))),'o',ibindmut(:,1),log(ibindmut(:,2)/max(ibindmut(:,2))),'x');shg
%figure(7);plot(ibind(:,1),ibind(:,2)/58,'o',ibindmut(:,1),ibindmut(:,2)/69,'x');shg