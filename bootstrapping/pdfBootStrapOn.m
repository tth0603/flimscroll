function pdfBootStrapOn(vector,numberOfBootStraps,bins,figNum)
%this will make a PDF with error bars that correspond to the std for each bin
%(symmetrical error bars). As I recall this script is balls slow, so start
%by choosing a small number of boot straps to test the script.
%
%INPUTS:
%   VECTOR == a column vector of the data in the PDF
%   numberOfBootStraps ==  a scaler of the number of desired boot strapped 
%                          data sets (start small)
%   bins == a vector of bin sizes (see USAGE for example format)
%   figNum == a scaler specifying the figure to plot in
%OUTPUTS:
%   It'll make a 
%USAGE:
%pdfBootStrapOn(vector,100,[0,25,50:50:150,200:100:400,600],50)

%add err bars to data pdf
events=vector;
nBS = numberOfBootStraps;
out = bootstrapFromVector(events,nBS);
n = length(events);
%bins = [0,25,50:50:150,200:100:400,600]; 
db=diff(bins);

dat=vector;  
n=length(dat);     
%bins = [0,25,50:50:150,200:100:400,600]; 
db=diff(bins);
counts=histc(dat,bins)/n;
normFactor=1./([db db(end)]).*counts';  %you may need to 'unprime' the counts vector
figure(figNum);hold on;bar(bins,normFactor,'histc');shg

bsMat = [];
for i = 1:nBS
    dat = out.bootstraps(:,i);
    counts=histc(dat,bins)/n;
    normFactorTwo=1./([db db(end)]).*counts';
    bsMat = [bsMat normFactorTwo'];
end 
L = length(normFactor);
for i = 1:L
    avgVect(i) = mean(bsMat(i,:));
end

for i = 1:L
    %stdVect(i) = std(bsMat(i,:))/sqrt(L); %this is SEM, uncomment this and comment out the next line to use 
    stdVect(i) = std(bsMat(i,:));   %this uses the std deviation
end
    
x=bins+[db bins(end)]*0.5;
figure(figNum);errorbar(x,normFactor,stdVect,'b.');



