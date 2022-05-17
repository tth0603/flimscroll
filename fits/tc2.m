function pc = tc2(params,data1,data2)
% Will test out a maximum likelihood algorithm for fitting a distribution
% function.  in this instance we will try just a simple single state stat mech probability with 
% fit params
%       r == trxn rate 
%       k1 == asso'n const 1 
%       k2 == asso'n const 2
% 
% rateOfmRNA1 = r*k1/(1+k1)
% rateOfmRNA2 = r*k2/(1+k2)
%       
%   INPUTS: 
%       params == a n = 3 vector of initial guess for r and k1, k2. Ex: inargs =
%       [5 10 20];
%       data1 == a vector of measured mRNA rates for condition 1
%       data2 == a vector of measured mRNA rates for condition 2
%
% call via:
%       
%       fminsearch('tc1',params,[],data1,data2)
%
if (sum(params < 0) > 1) % check if there is any negative number in input variable
    pc = 1e9;    % give a big value to the result
    return;                % return to fminsearch - do not execute the rest of the code
end %this doesn' twork for reasons that are not clear to me!

r = params(1);
k1 = params(2);
k2 = params(3);
guess1 = r*k1/(1+k1);
guess2 = r*k2/(1+k2);

probvect = abs(data1-guess1+data2-guess2); %how best to combine and compare these?
% probvect = abs(data1*guess1+data2*guess2); %bad
% probvect = abs(data1*guess1+data2*guess2); 
prodprob = sum(log(probvect)); 


pc = prodprob; %should this be negative? it depends on how you combine and compare above


%from LJF:
                                      % Take product of all probabilities;
                                      % The fminsearch will then minimize the
                                      % -prodprob (i.e. maximize prodprob)
                                      % so we will maximize the likelihood
                                      % that the intervals vector
                                      % represents the said distribution
                                      % (The sum of the log of the vector entries
                                      % is the same as the log of the
                                      % product of the vector entries.  We
                                      % use the sum( log() ) because taking
                                      % the product of the entries yields a
                                      % number too small for the computer
                                      % to handle.
