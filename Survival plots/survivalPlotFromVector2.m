    function [mib,dataHandle,fitHandle] = survivalPlotFromVector2(data,col,fit,tm,tx,figNum,varargin)
% [mib,plotHandle] = freqSurvivalPlotFromVector2(data,col,fit,tm,tx,figNum,varargin)
%
%   fn to plot a survival curve and optionally fit to a single or double
%   exponential.
%       calls fns expfallone_mxl & expfalltwo_mxl
%
%   v2 is same as survivalPlotFromVector, but you can specify the fit curve
%   color via fitHandle output
%
%   INPUTS
%   data = matrix containing data to be plotted
%   col = column within matrix you wish to plot cureve for. this data
%   should be time intervals of events
%   fit = '0' for no fit, '1' for single exp, '2' for double exp
%   tm = the time resolution of your data in s
%   tx = how long you ran the experiment in s
%   figNum = figure to plot in. Leave empty ('[]') if no plot, only a fit, is desired
%   varargin = the initial guesses for your fit params. only applies for
%   single and dbl expl fits. A single scaler for single expl (usu < 1, as
%   this is a guess of the rate, k[s^-1],). a vector length 3 for double 
%   exponential: [fraction in fast population(btwn 0 & 1), guess for fast
%   rate (less than one), guess for slow population rate (less than fast
%   guess].
%   
%   OUTPUT
%   mib == struc with fit rate constants and other suff and a plot of your data
%   against the fit
%   plotHandle == the handle for the surv. curve plot for 
%   specification of curve color, weight, etc using
%   set(plotHandle,'Color','k'), etc...
%
%   USAGE:
%       [struc,~] = survivalPlotFromVector(mat,1,2,4,3600,30,[0.8 0.1 0.01])
%       [~,p] = survivalPlotFromVector(ciaCell{i},7,0,30,3600,fignum); 
%       set(p,'Color',[0 0 0]);
%
%Harden 2014, updated 2018

% figNum = 29;

intervals=data(:,col);        %here we use time, not frames.
l=max(intervals);

avg=mean(intervals);
sigma=std(intervals);
n=length(intervals);

%sints=zeros([l,2]);        %allocating space, a vector of length l
sints=[];
for i=0:max(l)
    logik=intervals>i;
    sints=[sints;i sum(logik)];
end

if fit == 1
    %For single exp fit:
    tau=varargin{1};    %tau refers to rate constants, change to time constants by    
                 %edit expfalltwo_mxl       
    params=fminsearch('expfallone_mxl',tau,[],intervals,tm,tx);
    %corr=expfallone_mxl_correction(tau,intervals,tm,tx);
    rateConst=abs(params(1));  %gives rate constants for fit
    k1=1./rateConst;

    tfit=0:max(intervals);
    SPfit=exp(-tfit*rateConst);

    corr.Offset=0;      %if you decide to use Friedman's correction function comment these out
    corr.Nzero=max(sints(:,2));
    if ~isempty(figNum)
        figure(figNum);
        %plot(sints(:,1),log((sints(:,2)+corr.Offset)/corr.Nzero),'b',tfit,log(SPfit),'c');shg
        dataHandle = plot(sints(:,1),log((sints(:,2)+corr.Offset)/corr.Nzero),'LineWidth',1);
        hold on;
        fitHandle = plot(tfit,log(SPfit),'Color',[0.4 0.4 0.4]);shg
    else 
        dataHandle = 1; %do this only to avoid an error
    end
%     plot(sints(:,1),((sints(:,2)+corr.Offset)/corr.Nzero)); %,'b',tfit,log(SPfit),'c');shg
    mib.tConst=k1;
    mib.k1=rateConst;

%for double exponential fit:
elseif fit==2       %an example of a proper way to fit survival curves with a double exponential max likelihood
    inarg = varargin{1};
    
    params=fminsearch('expfalltwo_mxl',inarg,[],intervals,tm,tx); 
%     params=fminsearch('expfalltwo_mxl',inarg,[],sints(:,2),tm,tx); 
    %use this instead to account for a finite observation window:
    %params=fminsearch('TwoStepConvolution_mxl',inarg,[],intervals,tm,tx);   %calcs the rates and amplitudes of the mxl fit
%%%To correct for a finite observation time, uncomment 'corr' and comment
%%%'corr.Nzero' and 'corr.offset'%%%
    %corr=expfalltwo_mxl_correction(inarg,intervals,tm,tx);      %correction factors for the fact that want to characterize the entire phenomena, but only observe it for a part of the time
        %'corr' is a struc: 
        % corr.CorrectionFactor       % = (  1/( a*A + (1-a)*B ) 
        % corr.Nzero                  % = actual # of events in data set (observed+not observed)  
        % corr.Offset            % # remaining events at time = tx
        
    ap=params(1);
    a=1/(1+ap^2);
%     a=params(1);
    k1=params(2);
    k2=params(3);  %defines rate constants for fit

    tfit=[0:max(intervals)];

    SPfit=a*exp(-tfit*k1)+(1-a)*exp(-tfit*k2);  %This is the survival 
                                                %curve model
%     SPfit=k1*k2/(k1-k2)*(exp(-tfit*k1)-exp(-tfit*k2)); 

%this is the expn that we actually use to fit to account for exp
%resolution:
%     A=( exp(-tm*k1) - exp(-tx*k1) );
%     B = ( exp(-tm*k2) - exp(-tx*k2) );
%     SPfit=(  1/( a*A + (1-a)*B  ) )*...
%                                  ( a*k1*exp(-tfit*k1)+(1-a)*k2*exp(-tfit*k2) )+...
%                                  0;
%     SPfit = a*k1*exp(-tfit*k1)+(1-a)*k2*exp(-tfit*k2);

    corr.Offset=0;
    corr.Nzero=max(sints(:,2));
    if ~isempty(figNum)
        figure(figNum);
        dataHandle = plot(sints(:,1),log((sints(:,2)+corr.Offset)/corr.Nzero),'LineWidth',1);
        hold on;
        fitHandle = plot(tfit,log(SPfit),'Color',[0.4 0.4 0.4]);shg
        %plot(tfit,SPfit,'g-');
    else 
        dataHandle = 1; %do this only to avoid an error
    end
    mib.k1=params(2);
    mib.k2=params(3);
    mib.a = a;

%For no fit:
else
%     if xmax
    figure(figNum);
%      plot(sints(:,1),log((sints(:,2))/n));     %lin v log
%     plotHandle = stairs(sints(:,1),log((sints(:,2))/n),'LineWidth',1);
    dataHandle = stairs(sints(:,1),log((sints(:,2))/max(sints(:,2))),'LineWidth',1);
%     plotHandle = stairs(sints(:,1),sints(:,2),'LineWidth',1);
%     else
%         plot(sints(1:xmax,1),((sints(1:xmax,2))/n));
%     end

end

mib.avg=avg;
mib.sigma=sigma/sqrt(n);
mib.n=n;
mib.yints = sints(:,2)/n;
mib.xints = sints(:,1);
