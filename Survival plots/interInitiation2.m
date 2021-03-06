function output = interInitiation2(cia,lowLengthThrd,figureNumber)
%Use this function when simultaneously removing short high events and
%plotting interinitiation times.  Plug in the unedited cia, this returns a
%cia of purely 1 evnts.
%output = interInitiation2(cia,lowLengthThrd,figureNumber)
%this function plots the time between the start of successive 'hi' events.
%   cia: load CumulativeIntervalArray file
%   figureNumber: The fig number you are plotting in.
%   lowLengthThrd: threshold below which events will be excluded
%   output == (a structure with the following: 
%   [rateConstant average standardDeviation(sigma) numberOfEvents(n) intervals(a struc of the plotted event times) sints(a struc of the sorted event times)]
%   
logik=(cia(:,1)==1);        %find the up events
cia2=cia(logik,:);

logik2=(cia2(:,5)>lowLengthThrd);
cia3=cia2(logik2,:);        %get rid of the short events, a cia-like array
    
L=length(cia3);
ints=[];
for i=1:L-1     %make a vector of the times between successive initiations on the same AOI
    if cia3(i+1,7)==cia3(i,7)
        ints=[ints;cia3(i+1,2)-cia3(i,2)];
    else
        ints=ints;
    end
end

output.avg=mean(ints);
output.sigma=std(ints);
output.n=length(ints);

%%To fit with a two exponential:
%     ap=1;      %%%%%%%%%%%%%%%You will need to fuck wiht these numbers%%%%%%%%%%%%%%%%%%
%     tau1=0.5;
%     tau2=0.1;
% 
%     tm=10;      %%%%%%%%%%%%%%%You will need to fuck wiht these numbers%%%%%%%%%%%%%%%%%%
%     tx=max(ints);
%     inarg=[ap tau1 tau2];
%     
%     params=fminsearch('expfalltwo_mxl',inarg,[],ints,tm,tx);   %calcs the rates and amplitudes of the mxl fit
% %%%To correct for a finite observation time, uncomment 'corr' and comment
% %%%'corr.Nzero' and 'corr.offset'%%%
%     %corr=expfalltwo_mxl_correction(inarg,intervals,tm,tx);      %correction factors for the fact that want to characterize the entire phenomena, but only observe it for a part of the time
%         %'corr' is a struc: 
%         % corr.CorrectionFactor       % = (  1/( a*A + (1-a)*B ) 
%         % corr.Nzero                  % = actual # of events in data set (observed+not observed)  
%         % corr.Offset            % # remaining events at time = tx
%         
% ap=params(1);
% a=1/(1+ap^2);
% k1=params(2);
% k2=params(3);  %defines rate constants for fit
% 
% tfit=[0:max(ints)];
% 
% SPfit=a*exp(-tfit*k1)+(1-a)*exp(-tfit*k2);  %This is the survival 
%                                             %curve model
% m=max(ints);
% sints=[];           %Sort the interinitiaion times for a survival plot
% for i=1:m;
%     logik=ints>i;
%     sints=[sints;i sum(logik)];
% end

      %%%%%%%%%%%%%%%You will need to fuck wiht these numbers%%%%%%%%%%%%%%%%%%
tau=0.01;    %tau refers to rate constants, change to time constants by    
             %edit expfalltwo_mxl       
tm=10;      %%%%%%%%%%%%%%%You will need to fuck wiht these numbers%%%%%%%%%%%%%%%%%%
tx=max(ints);      %tm is about the frame length, of the lower event threshold 
                        %tx is more of less the longest possible observable
                        %event.

params=fminsearch('expfallone_mxl',tau,[],ints,tm,tx);
%corr=expfallone_mxl_correction(tau,intervals,tm,tx);

rateConst=params(1);  %gives rate constants for fit
timeConst=1./rateConst;

tfit=0:max(ints);
SPfit=exp(-tfit*rateConst);

sints=[];
for i=1:max(ints);
logik=ints>i;
sints=[sints;i sum(logik)];
end

corr.Offset=0;      %if you decide to use Friedman's correction function comment these out
corr.Nzero=max(sints(:,2));

figure(figureNumber);
plot(sints(:,1),log((sints(:,2)+corr.Offset)/corr.Nzero),'b',tfit,log(SPfit),'c');shg
% output.amplitude=a;
% output.k1=params(2);
% output.k2=params(3);
output.cia3=cia3;       %a cia-like array of up events
output.sints=sints;
end