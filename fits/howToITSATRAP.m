%script instructing on how to fit data with a trapazoid with "base" points
%'a' and 'd' and 'shoulder' points 'b' and 'c'. also will give you an
%amplitude 'amp' and background level 'bkg'. input guesses for these to fitsatrap and
%it will spit these fitted terms back at you in the form of a vector. then
%plug this into itsatrap and use that output to plot the fit curve against
%your data 

%make some test data (or use real data):
xmax = 100;
inarg = [10 50 15]; %amp mean std
r = rand(xmax,1); %noise term
ys = inarg(1)*exp( -( ([1:xmax]'-inarg(2)).^2/(2*inarg(3)^2) ) )+r; %gaussian dist with noise
%plot the test data:
figure(94);plot([1:xmax],ys);shg

%% the business:
xValues = [1:100]'; %same size as your data (in this case 'ys'), these are the x values of your data

%fit the data:
fitTerms = lsqcurvefit('itsatrap',[20 45 55 90 10 1],[],xValues,ys); %ys are your data. 'fitTerms' is the vector to plug into the trap plotting function [a b c d amp bkg]

%generate the fitted trap:
trap = itsatrap(fitTerms,xValues); %if your data is low res you may want to use a different set of xValues, like [1:0.01:100]

%plot everything:
figure(99);plot(xValues,trap);shg  %trp
hold on;plot([1:xmax],ys);shg  %the data that were fit


