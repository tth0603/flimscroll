function mih = ciaChopper(varargin)
% out = ciaChopper(analData1,analData2,experimentTime)
% 
% This script will omit all events which begin after a certain time,
% allowing for the computaion of event frequencies and comparison between
% experiements. That certain time should be the the shortest experiment
% (relative to the start of NC14). 

experimentTime = varargin{3};

% Exclude all events from each cia that begin after experimentTime:
ciaCell = {};
for i = 1:2
    cia = varargin{i}.cia;
    logi = cia(:,7) < experimentTime;
    ciaCell{i} = cia(logi,:);
end

%ficure out nuke numbers
nukeNumV = zeros(2,1);
for i = 1:2
    nukeNumV(i) = max(varargin{i}.nukesMat(:,3)); %third column is nuke column
end
%add spotsMat, nukesMat, and cia info from each subsequent data set to the
%first data set
for i = 2  
    %combine cia
    mat3 = ciaCell{i};
    mat3(:,1) = mat3(:,1) + sum(nukeNumV(1:i - 1));
    mih = [ciaCell{1}; mat3];
end
