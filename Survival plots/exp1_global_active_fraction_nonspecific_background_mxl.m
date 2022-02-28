function pc=exp1_global_active_fraction_nonspecific_background_mxl(inarg,datacells,rc)
%
%function exp1_global_active_fraction_nonspecific_background_mxl(inarg,datacells,rc)
%
% Intended to fit multiple data sets (of e.g. RNAP first landings on dsDNA)
% to single exponentials using a global active fraction.  In each instance
% it will account for the number of AOIs for which we do not observe any
% landing, and for the number of AOIs already occupied at time zero.
%
% inarg == [ap r1 r2 r3 ...rN]  where 
%         (active fraction) = ap^2/(1+ap^2)   and
%     r1, r2 r3 are the rates for the rising exponential fits of the form:
%  (active fraction)*(1 - exp(-t*ri) ) - afc*(1 - exp(-t*rc) )   (2nd term is nonspecific 
%                                               backgrgound surface binding)
% datacells == cell array, where each member contains one experimental data set  
%     i.e. i = 1 to N for N data sets where
% datacells{i}.intervals==vector list of first binding intervals measured for
%                   those AOIS in which at least one landing event was observed
% datacells{i}.tx == maximum observation time 
% datacells{i}.Nt == total number of AOIs, including those already
%                 occupied at t = 0 and those for which no landing was observed
% datacells{i}.Ns == number of AOIs for which no landing was observed
% datacell{i}.Nz == number of AOIs already occupied at time = 0
% rc == rate for nonspecific surface binding 
% call via fminsearch('exp1_global_active_fraction_nonspecific_background_mxl',inarg,[],datacells,rc)
%`z
% see notes from 5/25/2012 for derivations
% Also see notes from 6/2/2012
Af=inarg(1)^2/(1+(inarg(1))^2);     % Active fraction = Af  
DataNum=length(datacells);          % Number of data sets we will fit
prodprob=0;                         % Initialize likelihood function value
for indx=1:DataNum
    intervals= datacells{indx}.intervals;   % List of landings
    tx=datacells{indx}.tx;                  % max observation time
    Nt=datacells{indx}.Nt;                  % total AOI #\
    Ns=datacells{indx}.Ns;                  % AOI# w/o landing
    Nz=datacells{indx}.Nz;                  % AOI# occupied at t=0;
    rate=inarg(1+indx);                     % landing rate for this data set
                                    % Cycle through all the data sets
            % List of probabilities for the observed landings
    probints = ( (Af*Nt - Nz)/Nt )*(rate+rc)*exp(-(rate+rc)*intervals) + (1-Af)*rc*exp(-rc*intervals);   
        %      (#active past t=0)/(total AOI#) *rate*exp( )

        % Probability for observing one instance of AOI w/o a landing 
    probNs = ( ((Nt - Af*Nt)/Nt)*exp(-tx*rc) + ( (Af*Nt - Nz)/Nt )*exp(-tx*(rate+rc)) );   
            % (inactive #)/(total #)+  (# active w/o landing prior to tx)/(total AOI #)
       % (binds w/ nonspecific rate rc)    (  binds with rate of (specific rate)+(nonspecific rate rc)  ) 
         % Collect running log of product of the above probabilities.
         % Note the Ns factor to account for the Ns instances of AOIs 
         % that  are inactive wrt specific binding and instead bind only at
         % the nonspecific rate = rc (see above formulas)
    prodprob=prodprob + sum(log(probints)) + Ns*log(probNs);
end
                                                            

% maximize product of all probabilities;
pc=-prodprob;   
% 
% cumulative plot fit should look like:
%  (Af*Nt-Nz)*( 1 - exp(-t*(rate+rc)) ) +  Nz  +  Nt*(1-Af)*( 1-exp(-t*rc) )
%  (probably also divide above quantity by Nt for plot to approach 1)

			% Now we need to fabricate a data set with 
            % active fraction = 1 to force our fit of the 
			% control data to have an af = 1 as well
%%f=0.001:.001:.98;	% where f = ( 1=exp(-t/?) ), run out to a fraction 0.98
%%ints=-4066*log(1-f);	% Solve for the event times with tau = 4066
%%datacells{2}.intervals=ints;	% There are 980 interval times here.  Note that 
             %  980/.98=1000, so we have 20 fake AOIs not
             % yet occupied at our final time of max(ints) = 15906 s
%%datacells{2}.Nz=0;		% None of our fake AOIs occupied at t=0 
%%datacells{2}.tx=15906;
%%datacells{2}.Ns=20;		% 20  AOIs w/ no landing, as noted immediately above
%%datacells{2}.Nt=1000;	% 1000 total AOIs