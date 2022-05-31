%This will give survival curves for both the WT and the mutant landings.
%must have interval data file for each.  Start by loading the file:

load t:\matlab\data\b1p63a.dat -mat   %Creates 'Intervals' variable.

cia=Intervals.CumulativeIntervalArray;
logik=cia(:,1)==1;
ints=cia(logik,4);                  %These are the landing intervals

disp('Sample WT intervals:');
disp(ints(1:10));

sints=[];                        %creates the variable that will be the 
                                %sorted landing intervals
for indx=0:620
    logik=ints>indx;
    sints=[sints;indx sum(logik)];
end                             %Creats a vector of the sorted intervals
                                %Differs from the sorted intervals that Ben
                                %showed you in b1p64a.m.  This vector gives
                                %at each time interval (every second) gives
                                %the number of landings that are of that
                                %duration or longer
%Repeat for Mutant Data:
load t:\matlab\data\b1p62c.dat -mat

ciamut=Intervals.CumulativeIntervalArray;
logik=ciamut(:,1)==1;
intsmut=ciamut(logik,4);

disp('Sample MUT intervals:');
disp(intsmut(1:10));

sintsmut=[];

for indx=0:690
    logik=intsmut>indx;
    sintsmut=[sintsmut;indx sum(logik)];
end  

figure(5);plot(sints(:,1),log(sints(:,2)/221),'o',sintsmut(:,1),log(sintsmut(:,2)/290),'x');shg