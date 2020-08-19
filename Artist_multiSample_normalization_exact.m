function Artist_multiSample_normalization_exact(base,sampleID)
% last edit July-27-2020
%>> Artist_mutliSample_normalization_exact('lowest',{'mucin1','mucin2'})

% Changes from v1:
    % remove rounding & collect unrounded mean from controlsims
        % output from Artist_myRun
    % scaling after psuedocount
    % sample name specification and retrieval explicit in user input
    % output filename holds input sampleIDs
    
% Acts on 'MWUdone*' .mat files in currrent directory
% INPUTS
    % User-defined <base> for cross-sample normalization of read depth
        % 'cpm' or 'cp10m' == counts per 1 or 10 million
        % 'lowest' == sample with lowest total reads
    % <sampleID> is cell array of strings UNIQUE to sample filenames
% PROCESS
    % Takes mean of Artist normalized Inputs for each Input--Sample pair
    % Rounds to nearest integer
    % Adds pseudocount= 1 output count (unit = cpm, cp10m, or lowest)
% OUTPUT
    % TA-sites
    % mean of Artist normalized Inputs, also normalized by Sample depths
    % depth-normalized Samples
    
%% Load Artist output files in current directory
%input list of file names, held in structure
files = dir('MWUdone*'); 
%reorder from UNIX to natural numeric sequence
files = natsortfiles({files.name});
idx=contains(files,sampleID);
files=files(idx);

%% Consolidate sample names & read data into structure
s = struct('sampleName',{},'inData',[],'outData',[]);
depth=[];
for i=1:length(files)
    id = contains(files,sampleID(i));
    infile = char(files(id));
    s(i).sampleName = sampleID(i);
    load(infile,'myNormInput','mySample','myTAsites');
    s(i).inData = [myNormInput, mySample];
    %Add pseudocount to all
    s(i).outData = s(i).inData +1;
    %capture read depth for each Input-Sample pair
    depth(i)=sum(s(i).outData(:,2));
end

%% Calculations
% Normalize across samples to CPM (or CPM10M)
if strcmp(base,'cpm')==1
    scale=10^6;
elseif strcmp(base,'cp10m')==1
    scale=10^7;
% Normalize across samples to lowest read depth for sample set
elseif strcmp(base,'lowest')==1
    scale=min(depth); 
end
%Add pseudocount, applying scaling factor
inputs=[];samples=[];names={};
for i=1:length(depth)
    %s(i).outData = s(i).inData +1;
    s(i).outData = s(i).outData .* scale/depth(i);
    inputs(:,i)=s(i).outData(:,1);
    samples(:,i)=s(i).outData(:,2);
end
%Mean of inputs
inputMean=mean(inputs,2);

%% Export table
TAsites = num2cell(myTAsites);
inputMean = num2cell(inputMean);
samples = num2cell(samples);
table_out = [TAsites, inputMean, samples];

%dirName = pwd; 
%ind = strfind(dirName,'/');
%dirName = dirName(ind(end)+1:end);

Tout = cell2table(table_out,'VariableNames',[{'TAsites'},{'normInput'}, sampleID]);
outName = strcat(strcat(cell2mat(sampleID)),'_sampleDepthNormalized-',base);
writetable(Tout,outName,'FileType','text','Delimiter','tab');



