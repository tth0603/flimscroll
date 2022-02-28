function mib = combineAnalData(varargin)
%
% To combine multiple data sets made through flimscroll2. Changes of note
% are: 1. 'relative frame number' (col 11 in spots mat and nukes mat)
% becomes nonsensical. 2. we change the nucleus numbers of each nucleus in
% all data sets after the first.
%
% INPUTS:
%   any number of analData structures from flimscroll2
%   EX: combineAnalData(analData1,analData2);
% 
% OUTPUTS:
%   a single analData structure with cia, spotsMat, & nukesMat
%
% Harden 2020

%get number of data sets
datnum = length(varargin);
%initiate the output struc by throwing in all the data from the first
mib = varargin{1};
%ficure out nuke numbers
nukeNumV = zeros(datnum,1);
for i = 1:datnum
    nukeNumV(i) = max(varargin{i}.nukesMat(:,3)); %third column is nuke column
end
%add spotsMat, nukesMat, and cia info from each subsequent data set to the
%first data set
for i = 2:datnum  %skip the first data set as we already added that
    %combine spotsMat
    mat = varargin{i}.spotsMat;
    mat(:,3) = mat(:,3) + sum(nukeNumV(1:i - 1)); %to each nuke num, add the sum of the nukes from all previous data sets
    mib.spotsMat = [mib.spotsMat; mat]; %tack the spots mat of the ith data set onto the spots mat that we will output
    %combine nukes mat
    mat2 = varargin{i}.nukesMat;
    mat2(:,3) = mat2(:,3) + sum(nukeNumV(1:i - 1)); %to each nuke num, add the sum of the nukes from all previous data sets
    mib.nukesMat = [mib.nukesMat; mat2];
    %combine cia
    mat3 = varargin{i}.cia;
    mat3(:,1) = mat3(:,1) + sum(nukeNumV(1:i - 1));
    mib.cia = [mib.cia; mat3];
end
