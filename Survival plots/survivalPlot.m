function output = survivalPlot(cia,eventType,fitType,figNum)
%Plot a survival curve of any specified single event type. This, presently
%will not fit the curve with an expenential.
%   useage: ouput = survivalPlot(cia,eventType,figNum);
%   Inputs:
%       cia == cia like array
%       eventType == hi or lo events in the form of a scalar. exp: [1]
%       figNum == figure to plot to
%       fitType == specifies fit. '0' for no fit, '1' for single exp, '2' for
%       two exp
%   output struc:
%       .intervals == un-sorted list of the plotted ints
%       .avg
%       .sigma     == average and standard deviation of events
%       .n         == number of events
%   USAGE:
%   output = survivalPlot(test,[-3],0,30);
if length(eventType)==2
    logik1=(cia(:,1)==eventType(1)|cia(:,1)==eventType(2));
    data=cia(logik1,:);
elseif length(eventType)==3
    logik1=(cia(:,1)==eventType(1)|cia(:,1)==eventType(2)|cia(:,1)==eventType(3));
    data=cia(logik1,:);
else
    logik1=cia(:,1)==eventType;
    data=cia(logik1,:);
end
intervals=data(:,5);        %here we use time, not frames.
l=max(intervals);

avg=mean(intervals);
sigma=std(intervals);
n=length(intervals);

%sints=zeros([l,2]);        %allocating space, a vector of length l
sints=[];
for i=1:max(l);
logik=intervals>i;
sints=[sints;i sum(logik)];
end

%%%%%%%for no fit%%%%%%%%%%
if fitType == 0;
figure(figNum);
%plot(sints(:,1),(sints(:,2)),'b');shg
plot(sints(:,1),log((sints(:,2))/n),'b');shg

%%%Fit with single exponential%%%%
elseif fitType == 1;
tau=0.001;    %tau refers to rate constants, change to time constants by    
             %edit expfalltwo_mxl       
tm=1;      %%%%%%%%%%%%%%%You will need to fuck wiht these numbers%%%%%%%%%%%%%%%%%%
tx=max(intervals);      %tm is about the frame length, of the lower event threshold 
                        %tx is more of less the longest possible observable
                        %event.
params=fminsearch('expfallone_mxl',tau,[],intervals,tm,tx);
%corr=expfallone_mxl_correction(tau,intervals,tm,tx);

rateConst=params(1);  %gives rate constants for fit
k1=1./rateConst;

tfit=0:max(intervals);
SPfit=exp(-tfit*rateConst);

corr.Offset=0;      %if you decide to use Friedman's correction function comment these out
corr.Nzero=max(sints(:,2));

figure(figNum);
plot(sints(:,1),log((sints(:,2)+corr.Offset)/corr.Nzero),'b',tfit,log(SPfit),'c');shg
output.k1=k1;


%%%%To fit with a two exponential%%%%
else %(fitType == 2;)
    ap=1.2;      %%%%%%%%%%%%%%%You will need to fuck wiht these numbers%%%%%%%%%%%%%%%%%%
    tau1=0.01;
    tau2=0.01;

    tm=2;      %%%%%%%%%%%%%%%You will need to fuck wiht these numbers%%%%%%%%%%%%%%%%%%
    tx=max(intervals);
    inarg=[ap tau1 tau2];
    
    params=fminsearch('expfalltwo_mxl',inarg,[],intervals,tm,tx);   %calcs the rates and amplitudes of the mxl fit
%%%To correct for a finite observation time, uncomment 'corr' and comment
%%%'corr.Nzero' and 'corr.offset'%%%
    %corr=expfalltwo_mxl_correction(inarg,intervals,tm,tx);      %correction factors for the fact that want to characterize the entire phenomena, but only observe it for a part of the time
        %'corr' is a struc: 
        % corr.CorrectionFactor       % = (  1/( a*A + (1-a)*B ) 
        % corr.Nzero                  % = actual # of events in data set (observed+not observed)  
        % corr.Offset            % # remaining events at time = tx
 
corr.Offset=0;      %if you decide to use Friedman's correction function comment these out
corr.Nzero=max(sints(:,2));

ap=params(1);
a=1/(1+ap^2);
k1=params(2);
k2=params(3);  %defines rate constants for fit

tfit=[0:max(intervals)];

SPfit=a*exp(-tfit*(k1))+(1-a)*exp(-tfit*(k2));  %This is the survival 
                                            %curve model
%SPfit=0.6*exp(-tfit*(1/150))+0.4*exp(-tfit*(1/1100)); 
figure(figNum);
%plot(sints(:,1),log((sints(:,2))/n),'b',tfit,log((1/corr.CorrectionFactor)*SPfit),'c');shg
plot(sints(:,1),log((sints(:,2)+corr.Offset)/corr.Nzero),'b',tfit,log(SPfit),'c');shg
output.k1=k1;
output.k2=k2;
output.a=a;
end

output.avg=avg;
output.sigma=sigma/sqrt(n);
output.n=n;
output.intervals=intervals;
output.sints=sints;
end