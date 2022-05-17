function pc = fitsatrap(inarg,data) %,xValues)
%itsatrap TH's version of matlab's Trapezoidal membership function.
%   ITSATRAP(INARG,DATA) returns a matrix which is the trapezoidal
%   membership function evaluated at X. PARAMS = [A B C D AMP BKG] is a 6-element
%   vector that determines the break points of this membership function and the amplitude and background.
%   We require that A <= B and C <= D. If B >= C.
%
%
%   See also TRAPMF, DSIGMF, EVALMF, GAUSS2MF, GAUSSMF, GBELLMF, MF2MF, PIMF, PSIGMF,
%   SIGMF, SMF, TRIMF, ZMF.

%   Roger Jang, 6-28-93, 10-5-93, 4-14-94.
%   Copyright 1994-2002 The MathWorks, Inc. 
%   UPDATED TH 2018

if nargin ~= 2
    error('Two arguments are required by the trapezoidal MF.');
elseif length(inarg) < 6
    error('The trapezoidal MF needs at least six parameters: a, b, c, d, amp, bkg');
end

a = inarg(1); b = inarg(2); c = inarg(3); d = inarg(4);
amp = inarg(5);
bkg = inarg(6);

if a > b,
    error('Illegal parameter condition: a > b');
elseif c > d,
    error('Illegal parameter condition: c > d');
end

y1 = zeros(size(data));
y2 = zeros(size(data));

% Compute y1
index = find(data >= b);
if ~isempty(index),
    y1(index) = ones(size(index))*amp;
end
index = find(data < a);
if ~isempty(index),
    y1(index) = ones(size(index))*bkg;
end
index = find(a <= data & data < b);
if ~isempty(index) & a ~= b,
    y1(index) = (data(index)-a)*(amp-bkg)/(b-a)+bkg;
end

% Compute y2
index = find(data <= c);
if ~isempty(index),
    y2(index) = ones(size(index))*amp;
end
index = find(data > d);
if ~isempty(index),
    y2(index) = ones(size(index))*bkg;
end
index = find(c < data & data <= d);
if ~isempty(index) & c ~= d,
    y2(index) = (d-data(index))*(amp-bkg)/(d-c)+bkg;
end

% Compute y
pc = min(y1, y2);

%saving this:

% %
% % This function will be called from using the
% % lscurvefit() function to a trapazoid.  The form of the
% % function is weird, but the parameters are: a and d locate the “feet” of 
% % the trapezoid and the parameters b and c locate the “shoulders.”
% %
% % so a < b < c < d
% %
% % the input arguements are given in inarg according to
% % inarg  == [ a b c d ]
% % data == x values which are the centers of the bins for the histogram
% %
% %Call using:
% %   out =lsqcurvefit('itsatrap',[2 3 7 9],xValues,yValues);
% %       where yValues are the probabilities from teh histogram
% %   out == a vecter of a, b, c, & d
% xValues = [1:0.1:100]';
% 
% pc = trapmf(xValues,[inarg(1) inarg(2) inarg(3) inarg(4)]);   
% 
% %to fit data with a gaussian and plot:
% %out =lsqcurvefit('gaussFitFunc',[10 0.0001 0.1],pointXvalues,normFactor);
% % x = [-0.3:0.001:0.3];
% % fitCurve = out(1) * exp(-( (x-out(2)).^2/(2*out(3)^2) ) );
% % hold on; plot(x,fitCurve,'r-');shg
% 



