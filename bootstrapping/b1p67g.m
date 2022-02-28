%This creates bootstraps for the REACTION rate constants. It needs an
%input of bootstrapped OBSERVED rate constants.  Examples of the input are
%b1p67a and b1p67f scripts, corresponding to the MUT and WT data from the
%29 July 2011 Experiment.

%Additional imput needed:

%K1=K1=3.1e-3/0.2e-9       %This is the 1st order rate constant. (ie 2nd
                           %order rate const multiplied by [RNAP].
%km1=0.1;                  %These starting parameters seem to work OK, but
%k2=0.1;                   %you may need to tweak them a bit.
%km2=0.01;

argrates=zeros(10,3);       %allocates space for the output

for i=1:100                 %THis loop takes the amplitude parameter that the 
b(i)=1/(1+argouts(i,1)^2);  %observed rates bootstraps spit out, which is in Larry's weird format, and 
end
a=b';                       %rearranges so it is in the form to plug into the eqn.

for i=1:100                 
longtau=1/argouts(i,3);
shorttau=1/argouts(i,2);
short2longratio=a(i)/(1-a(i)^2);
argrates(i,:)=fsolve('sigma54_kinetics_fit_values_r12_weight_ratio',[km1,k2,km2],[],K1,longtau,shorttau,short2longratio);
end

%This algorithm may return a list of n errors:
%"Optimization terminated: first-order optimality is less than options.TolFun."
%However, the output is still in 'argrates' and seems to be reasonable.