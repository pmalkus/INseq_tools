function [bootstrapcontrol]=simulateEqualSaturation_v3(inputReads,sampleReads,boots)
% PNM, May 2020: 
    % break down #boots to avoid memory overload; loop to total
    % simplified input, do calculations in script, not function call

%taproportion is the gap in ta density between libraries
%reads is the number of reads to simulate.  Reads should equal the total
%number of reads in the experimental condition that you want to compare to
%input reads is the vector of control input
% calculations to simplify input (PNM)
taproportion=length(sampleReads(sampleReads~=0))/length(inputReads(inputReads~=0));
reads=sum(sampleReads);
inputreads=inputReads;
% and carry on with original...
inputreadssum=sum(inputreads);
inputproportion=inputreads./inputreadssum;
inputproportiontanorm=inputproportion.*taproportion;
inputproportiontanorm(length(inputproportiontanorm)+1,1)=1-taproportion;

% added loop (PNM)
if boots>10
    %build from loop of smaller sets of simulations
    i=0; j=10; multinominputsample = [];
    % <i> number of simulations run; <j> number of simulations to run this loop
    while i<boots
    temp = mnrnd(reads,inputproportiontanorm,j);
    multinominputsample = [multinominputsample; temp];
    i=i+j;
        if i+j>boots
            j=boots-i;
        end
    end
else multinominputsample=mnrnd(reads,inputproportiontanorm,boots);
end

multinominputsample=multinominputsample';
multisum=sum(multinominputsample);
difference=multisum-multinominputsample(length(multinominputsample),1);
correction=repmat((reads./difference),length(multinominputsample),1);
correctedinput=multinominputsample.*correction;
bootstrapcontrol=correctedinput(1:length(inputreads),:);
end



