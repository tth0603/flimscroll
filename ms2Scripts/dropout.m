function mib = dropout(cia)
% mib = dropout(cia)
% 
% This script will take in a cia and get rid of singe frame signal
% dropouts.
% 
% INPUT is a cia, cumulative interval array, from analData via flimscroll2.
% Its columns are:
% '1.nucleus 2.first,middle,or Last(= -3,1,or 3) 3.frameStartRelativeToNCstart 4.frameEnd 5.deltaFrames 6.deltaTime(sec) 7.startTimeRelativeToNCstart'
% 
% OUTPUT is a cia.
% 
% Timothy Harden 2020

%first, figure out the acquisition time. 
% Aside: This is a bit complicated as we
% (somewhat arbitrarily) give single frame events half acquisition time
% intervals in flimscroll2. If you have complaints about that think about it
% and propose a better solution. I bet there is one but it's not trivial. 
%Anyways, the length of a dropout will be the full acquisiton time.
%So here we'll find a single frame event, then double it to come up with
%acqusition time:
singles = cia(cia(:,5) == 1,6);
acqnT = mean(singles)*2;

%get the nukes with signal:
nukes = unique(cia(:,1));
mib = [];
for i = nukes'
    mat = cia(cia(:,1) == i,:); %all rows pertaining to nuke i
    L = size(mat,1);
    if L > 1  %pick only nukes that have more than one event
        for j = 1:L - 1
%             try %for TSing
            if mat(j,4) + 2 == mat(j + 1,3) %does event j end one frame prior to the start of event j + 1?
                %if so, shove all of event j into event j + 1
                mat(j + 1,2) = mat(j,2); %give j's binary label to event j + 1
                mat(j + 1,3) = mat(j,3); %j + 1's start time is now j's start time
                mat(j + 1,5) = mat(j,5) + mat(j + 1,5) + 1; %sum their delta frame times, plus one for the dropout
                mat(j + 1,6) = mat(j,6) + mat(j + 1,6) + acqnT; %%sum their delta times, plus the time of acquisition
                mat(j + 1,7) = mat(j,7); %make the j + 1's start time that of j
                mat(j,:) = NaN(1,1); %make event j NaN for the moment
            end
%             catch
%                 keyboard
%             end
        end
    end
    mib = [mib; mat];
end
%get rid of NaNs:
logi = ~isnan(mib(:,1));
mib = mib(logi,:);

                
    

