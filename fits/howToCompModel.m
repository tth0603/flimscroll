%this is just a explanatory script for how to use mxl methods to determin
%params in CS's shadow enhancer data when we apply the competition model.
%We start with a couple of test cases ('tci') for proof of principle
% this calls a number of custom written fns for each case, namely:
%       pc = tc1(params,data)
%       pc = tc2(params,data1,data2)

%% for tc1:
%gen randon, norm dist data
n = 100;
r = 10+2*randn([n 1]); %r = 10 with some sorm dist noise
k = 5+randn([n 1]);  %k = 5 with norm noise
%inspect these params:
figure(30);hist(r);shg
figure(31);hist(k);shg
%now use these to get simulated data:
ms = r.*k./(1+k);

%to evaluate:
[mxlParams,fval,exitflag] = fminsearch('tc1',[1 2],[],ms)  %guesses are r = 1, k = 2. mxlParams = [mxlr mxlk]
%or
options = optimset('MaxFunEvals',1e6);
[mxlParams,fval,exitflag] = fminsearch('tc1',[1 2],options,ms);

%neither of these will converge, but it gets the correct (ish) answer
%anyhow. any by correct answer, I do not mean the correct input params but
%some set of params such that they describe the simulated data just as well
%as the input params. whether or not we get the real input params does not
%depend on the size of the simulated data we input.

%test how good
an = mean(ms)
guess = mxlParams(1).*mxlParams(2)./(1+mxlParams(2))

%% tc2
%gen randon, norm dist data
n = 100;
%big values:
r = 5+2*randn([n 1]); %r = 10 with some sorm dist noise
k1 = 500+randn([n 1]);  %k1 = 5 with norm noise
k2 = 5+randn([n 1]);  %k2 = 20
%small values (<1):
r = 5+2*randn([n 1]); %r = 10 with some sorm dist noise
k1 = 0.5+0.01*randn([n 1]);  %k1 = 5 with norm noise
k2 = 0.01+0.005*randn([n 1]);  

%now use these to get simulated data:
ms1 = r.*k1./(1+k1);
ms2 = r.*k2./(1+k2);

%to evaluate:
%[mxlParams,fval,exitflag] = fminsearch('tc2',[1 2 3],[],ms1,ms2);  %guesses are r = 1, k = 2. mxlParams = [mxlr mxlk]
%or
options = optimset('MaxFunEvals',1e9,'TolX',0.001); %may want to change the tolerance. not sure what the default is, or what this does
%options = optimset('MaxFunEvals',1e9);
[mxlParams,fval,exitflag] = fminsearch('tc2',[1 0.1 0.01],options,ms1,ms2);

%in future download and try fminsearchbnd for bounded params

%neither of these will converge, but it gets the correct (ish) answer
%anyhow. any by correct answer, I do not mean the correct input params but
%some set of params such that they describe the simulated data just as well
%as the input params. whether or not we get the real input params does not
%depend on the size of the simulated data we input.

%test how good
an1 = mean(ms1)
guess1 = mxlParams(1).*mxlParams(2)./(1+mxlParams(2))
an2 = mean(ms2)
guess2 = mxlParams(1).*mxlParams(3)./(1+mxlParams(3))

%% tc2 with bound params
%gen randon, norm dist data
n = 100;

%small values (<1):
r = 0.1+0.01*randn([n 1]); %r = 10 with some sorm dist noise
k1 = 0.5+0.01*randn([n 1]);  %k1 = 5 with norm noise
k2 = 0.01+0.005*randn([n 1]);  

%now use these to get simulated data:
ms1 = r.*k1./(1+k1);
ms2 = r.*k2./(1+k2);

%to evaluate:
%[mxlParams,fval,exitflag] = fminsearch('tc2',[1 2 3],[],ms1,ms2);  %guesses are r = 1, k = 2. mxlParams = [mxlr mxlk]
%or
options = optimset('MaxFunEvals',1e9,'TolX',0.001); %may want to change the tolerance. not sure what the default is, or what this does
%options = optimset('MaxFunEvals',1e9);
[mxlParams,fval,exitflag] = fminsearchbnd('tc2',[1 0.1 0.01],[0 0 0],[inf inf inf],options,ms1,ms2);
%  usage: x=fminsearchbnd(fun,x0,LB,UB,options,p1,p2,...)
%from: https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd--fminsearchcon

%in future download and try fminsearchbnd for bounded params

%neither of these will converge, but it gets the correct (ish) answer
%anyhow. any by correct answer, I do not mean the correct input params but
%some set of params such that they describe the simulated data just as well
%as the input params. whether or not we get the real input params does not
%depend on the size of the simulated data we input.

%test how good
an1 = mean(ms1)
guess1 = mxlParams(1).*mxlParams(2)./(1+mxlParams(2))
an2 = mean(ms2)
guess2 = mxlParams(1).*mxlParams(3)./(1+mxlParams(3))

%try systematically varying the initial guesses:
v = 0.01:0.01:1;
len = length(v);
mxlParams = zeros(len,3);
for i = 1:len
    mxlParams(i,:) = fminsearchbnd('tc2',[i i i],[0 0 0],[10 10 10],options,ms1,ms2);
end
%eval the means from this:
meanParams = mean(mxlParams);
an1 = mean(ms1)
guess1 = meanParams(1).*meanParams(2)./(1+meanParams(2))
an2 = mean(ms2)
guess2 = meanParams(1).*meanParams(3)./(1+meanParams(3))

%%try systematically varying the initial guesses independently:

%% the real mccoy:
r1 = 1+0.1*randn([n 1]); %r = 10 with some sorm dist noise
r2 = 0.1+0.01*randn([n 1]); %r = 10 with some sorm dist noise
k1d = 1+0.1*randn([n 1]);
k1p = 0.5+0.05*randn([n 1]);
k2d = 0.01+0.001*randn([n 1]);
k2p = 0.1+0.01*randn([n 1]);

figure(30);hist(r2);shg

m1d = r1.*k1d./(1+k1d);
m2d = r2.*k2d./(1+k2d);
m1p = r1.*k2p./(1+k2p);
m2p = r2.*k2p./(1+k2p);
m12 = (r1.*k1d+r2.*k2p)./(1+k1d+k2p);
m21 = (r2.*k2d+r1.*k1p)./(1+k2d+k1p);

params = [1 1 0.1 0.1 0.1 0.1];

options = optimset('MaxFunEvals',1e9);
[mxlParams,fval,exitflag] = fminsearchbnd('rateRaiderv1',params,[0 0 0 0 0 0],[inf inf inf inf inf inf],options,m1d,m2d,m1p,m2p,m12,m21);
% compare mxlParams with:
% params == [r1 r2 k1d k1p k2d k2p]

%test how good
an1 = mean(m1d)
guess1 = mxlParams(1).*mxlParams(3)./(1+mxlParams(3))
an2 = mean(m12)
guess2 = (mxlParams(1).*mxlParams(3)+mxlParams(2).*mxlParams(6))./(1+mxlParams(3)+mxlParams(6))


%try to get the errors on the estimates:
%[x,fval,exitflag,output,grad,hessian] = fminunc(___)














