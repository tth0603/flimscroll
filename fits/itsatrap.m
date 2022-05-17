function pc = itsatrap(inarg,x) %,xValues)
%itsatrap TH's version of matlab's Trapezoidal membership function.
%   ITSATRAP(INARG,DATA) returns a matrix which is the trapezoidal
%   membership function evaluated at X. INARG = [A B C D AMP BKG] is a 6-element
%   vector that determines the break points of this membership function and the amplitude and background.
%   We require that A <= B and C <= D. If B >= C.
%
%
%   See also TRAPMF, DSIGMF, EVALMF, GAUSS2MF, GAUSSMF, GBELLMF, MF2MF, PIMF, PSIGMF,
%   SIGMF, SMF, TRIMF, ZMF.

%   Roger Jang, 6-28-93, 10-5-93, 4-14-94.
%   Copyright 1994-2002 The MathWorks, Inc. 
%   UPDATED Harden 2018

if nargin ~= 2
    error('Two arguments are required by the trapezoidal MF.');
elseif length(inarg) < 6
    error('The trapezoidal MF needs at least six parameters: a, b, c, d, amp, bkg');
end

a = inarg(1); b = inarg(2); c = inarg(3); d = inarg(4);
amp = inarg(5);
bkg = inarg(6);

if a > b
    error('Illegal parameter condition: a > b');
elseif c > d
    error('Illegal parameter condition: c > d');
end

y1 = zeros(size(x));
y2 = zeros(size(x));

% Compute y1
index = find(x >= b);
if ~isempty(index)
    y1(index) = ones(size(index))*amp;
end
index = find(x < a);
if ~isempty(index)
    y1(index) = ones(size(index))*bkg;
end
index = find(a <= x & x < b);
if ~isempty(index) && a ~= b
    y1(index) = (x(index)-a)*(amp-bkg)/(b-a)+bkg;
end

% Compute y2
index = find(x <= c);
if ~isempty(index)
    y2(index) = ones(size(index))*amp;
end
index = find(x > d);
if ~isempty(index)
    y2(index) = ones(size(index))*bkg;
end
index = find(c < x & x <= d);
if ~isempty(index) && c ~= d
    y2(index) = (d-x(index))*(amp-bkg)/(d-c)+bkg;
end

% Compute y
% pc = min(y1, y2);
pc = y1 + y2;





