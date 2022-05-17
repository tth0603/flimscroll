%This will give initial binding curves for both the WT and the mutant landings.
%must have interval data file for each.  Start by loading the interval file:

load t:\matlab\data\b1p63a.dat -mat   %Creates 'Intervals' variable.

cia=Intervals.CumulativeIntervalArray;
logik=(cia(:,1)==-2);  %This logic variable corresponds to ALL on rate intervals   
%logik=(cia(:,1)==-2)|(cia(:,1)==0);
ints=cia(logik,4);                  %These are the times before initial binding

ibind=[];                        %creates the variable that will be the 
                                %sorted initial binding times
for indx=0:2690
    logik=ints<=indx;
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
logik=(ciamut(:,1)==-2);
%logik=(ciamut(:,1)==-2)|(ciamut(:,1)==-0);
intsmut=ciamut(logik,4);

ibindmut=[];

for indx=0:2690
    logik=intsmut<=indx;
    ibindmut=[ibindmut;indx sum(logik)];
end  

%MODEL for WT:

tau=0.0022;    %tau refers to rate constants, change to time constants by    
             %edit expfalltwo_mxl       
tm=0.5;
tx=3600;      %tm is (1/2)*frame length, tm is the longest event 
                        %in the data set

params=fminsearch('expfallone_mxl',tau,[],ints,tm,tx);

k=params(1);  %gives rate constants for fit

disp('WT Rate Constant (1/s):');    %Displays rate and time constants
disp(k*(1/0.2e-9));

disp('WT Time Constant (1/s):');
disp(1./k);

tfit=0:0.5:2690;

SPfit=exp(-tfit*k);  %This is the survival curve model

%ERROR for WT rate constants:

tbl=[[ints] 1/82*ones(82,1)];
indx=probability_steps(tbl,82);
argouts=zeros(100,1);

for i=1:100
    ind=probability_steps(tbl,82);
    argouts(i)=fminsearch('expfallone_mxl',tau,[],ind,0.5,3400);
end

k=(1/0.2e-9)*argouts;       %changing 2nd O to 1st O rate constant using exp'l [RNAP]
disp('WT Rate Constant Error (1/sM))');
disp(std(k));

%WT botstrap curves:


%btstrp=[];

%for i=0:2690
    %logik=indx<=i;
   % btstrp=[btstrp;indx sum(logik)];
%end

%MODEL for MUT:

tau=0.0031;    %tau refers to rate constants, change to time constants by    
             %edit expfalltwo_mxl       
tm=0.5;
tx=3600;      %tm is (1/2)*frame length, tm is the longest event 
                        %in the data set

paramsmut=fminsearch('expfallone_mxl',tau,[],intsmut,tm,tx);

kmut=paramsmut(1);  %gives rate constants for fit

disp('MUT Rate Constant (1/s):');    %Displays rate and time constants
disp(kmut);

disp('MUT Time Constant (1/s):');
disp(1./kmut);

SPfitmut=exp(-tfit*kmut);  %This is the survival curve model

%ERROR for MUT rate constants:

tbl=[[intsmut] 1/80*ones(80,1)];
indx=probability_steps(tbl,80);
argouts=zeros(100,1);

for indx=1:100
    ind=probability_steps(tbl,80);
    argouts(indx)=fminsearch('expfallone_mxl',tau,[],ind,0.5,3400);
end
   
k=(1/0.2e-9)*argouts;
disp('MUT Rate Constant Error (1/sM):');
disp(std(k));
                               
%PLOTS:
%fraction of initial binding events v. time:
figure(6);plot(ibind(:,1),ibind(:,2)/length(ints),'bo',ibindmut(:,1),ibindmut(:,2)/length(intsmut),'rx');shg 
%ln(unbound AOIS) v. timt:
figure(9);plot(ibind(:,1),log(1-ibind(:,2)/length(ints)),'bo',tfit,log(SPfit),'-',ibindmut(:,1),log(1-ibindmut(:,2)/length(intsmut)),'rx',tfit,log(SPfitmut),'-');shg
%,sboots(:,1:100),SP-1

%,btstrp(:,1),log(1-btstrp(:,2:20)/length(bootsam)),'-'