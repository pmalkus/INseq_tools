function INseq_read_filter_v3(infile)
%Last edit: June-24-2019 (changed cutRat default to 10^9)

%Tool for filtering sequencing data from Akk::Tn plate pools
    %low complexity sample & high coverage sequencing

%(Derived from "INseq_coveragePlot_v4.m")
    %Change file naming conventions to match "pipeline_v3.pl"
    %Preserve Goodman naming and structure, with Matlab filter optional
    %Make new file to preserve original/unfiltered read table
    %Remove user-set minimum for L and R reads
        %note: L or R can't =0 for L/R ratio calculation
    
% Process primary output from INseq pipeline
    %Apply cutoffs, plot, revise
    %Number of unique insertions, and matrix of coordinate/reads
    %Output file can be used for second-set of mapping-jobs (_SAMPLE.job2)

% INFILE is read mapping stats from from INseq pipeline: 
  %"INSEQ_experiment.scarf_<SAMPLE>.bowtiemap_processed.txt_<GENOME>"
  % preserve as "INSEQ_experiment.scarf_<SAMPLE>.bowtiemap_processed.txt_<GENOME>_RAW"
    % Designed for analysis of single sample
        %DOES NOT HANDLE concatenated input files
    % first column is genome name
    % second column is insertion position
    % third, fourth and fifth col = L, R, and Tot reads mapped to insertion
% OUTFILE is same format, saved as infile name:
  %"INSEQ_experiment.scarf_<SAMPLE>.bowtiemap_processed.txt_<GENOME>"


%Protection against overwriting RAW with filtered
inraw = strcat(infile,'_RAW');
files = ls; %FIX for Windows
if contains(files,inraw)==1
    disp('SCRIPT TERMINATED!  Use RAW data file as input for filtering');
    return
end
%Save a copy of INFILE (original) data with a new name, _RAW 
    %and handle input of _RAW for refiltering
raw = strfind(infile,'_RAW');   
if isempty(raw)==1
    inraw = strcat(infile,'_RAW');
    copyfile(infile, inraw);
    passfile = infile;
else
    passfile = infile(1:raw-1);
end

%initial/default cutoffs
cutTot = 1; %total number of reads (keep >= cutTot)
cutRat = 10^9; %ratio of L/R reads (keep <= cutRat)
minLR = 1; %number of L or R reads (keep >= minLR)
           %has to be >0 for ratio calculation!

%Retrieve sample name for file output
sampleStart = strfind(infile,'scarf_') +6;
sampleEnd = strfind(infile,'.bowtie') -1;
sampleName = infile(sampleStart:sampleEnd);

%Temporarilty save original file with a new name understood by 'readtable'
ind1 = strfind(infile,'_experiment') +10;
tempName = strcat(infile(1:ind1),'.txt');
copyfile(infile, tempName);

%read numeric columns from infile, sort by insertion coordinate
clrt = dlmread(tempName,'','B1..E(end)');%change '' to '\t' for tab only
clrt = sortrows(clrt,1);
%read character column from infile and store as cell array for output
T = readtable(tempName, 'ReadVariableNames',false, 'Delimiter', 'tab');
chrID = table2array(T(:,1)); %resize before output
%NOTE: only works with a single genome in sample
delete(tempName);

%FILTERING
%remove rows with totRead < initial cutoff for total reads (=1)
carCut = clrt(clrt(:,4) >= cutTot, :);
%remove rows with L or R < default cutoff (+1)
carCut_noz = carCut((carCut(:,2)>=minLR),:);
carCut_noz = carCut_noz((carCut_noz(:,3)>=minLR),:);
%calculate ratio of L-to-R, apply default cutoff (=1000)
ratLR=carCut_noz(:,2)./carCut_noz(:,3);
ccnzr = [carCut_noz ratLR];
ccnzr = ccnzr((ccnzr(:,5)>=(1/cutRat)),:);
ccnzr = ccnzr((ccnzr(:,5)<=(cutRat)),:);

%Plot total reads vs. ratio of LtoR, using DEFAULT CUTOFFs
f1 = figure;
scatter(ccnzr(:,4),ccnzr(:,5));
set(gca,'Xscale','log','Yscale','log');
xl=xlim; yl=ylim; %use later to preserve scale
set(gca, 'FontSize', 12)
name=strcat(sampleName,' (cutTotRead= ',num2str(cutTot),...
    ', cutLRratio=',num2str(cutRat),')');
title(name);
xlabel('total mapped reads (per coordinate)');
ylabel('ratio of L-to-R mapped reads (per coordinate)');
figfile = strcat(sampleName,'_filterV3_tot',num2str(cutTot),...
    '_LRrat',num2str(cutRat),'.fig');


%%
%Query user to revise cutTot, cutRat, and minLR values
  %Apply, then output coverage map, #of insertions passing filter
  %Save variable () to file for normalization, annotation, etc...
while(1)
    display('current value of cutTot='); cutTot,
    cutTot = input('Enter new cutoff for # of Total mapped reads per coordinate  ');
    display('current value of cutRat='); cutRat,
    cutRat = input('Enter the cutoff for LR ratio (high only) per coordinate  ');
    %display('current value of minLR='); minLR,
    %minLR = input('Enter new cutoff for # of LorR mapped reads per coordinate  ');
%create new matrix ('CTcut') to preserve info in 'ccnzr'
    CTcut = ccnzr((ccnzr(:,5) >= (1/cutRat)),:); 
    CTcut = CTcut((CTcut(:,5) <= cutRat),:);  
    CTcut = CTcut((CTcut(:,4) >= cutTot), :);
    CTcut = CTcut((CTcut(:,2) >= minLR), :);
    CTcut = CTcut((CTcut(:,3) >= minLR), :);
    aboveCut=size(CTcut,1);
    msg=strcat('Number of insertion sites above cutoffs (totReads>',num2str(cutTot),...
        ', LRratio<',num2str(cutRat),')  =',num2str(aboveCut));
    display(msg);
    %display('Approximate number of unique clones'), num2str(aboveCut);

    %plot total reads vs. ratio of LtoR, AFTER FILTERING
    f2 = figure;
    scatter(CTcut(:,4),CTcut(:,5)); %
    set(gca,'Xscale','log','Yscale','log');
    xlim(xl); ylim(yl); %apply scale from first plot to subsequent plots
    set(gca, 'FontSize', 12);
    name=strcat(sampleName,' (FILTERS:cutTot=',num2str(cutTot),...
        ', cutRat=',num2str(cutRat),')');
    title(name);
    xlabel('total mapped reads (per coordinate)','FontSize',14);
    ylabel('ratio of LtoR mapped reads (per coordinate)','FontSize',14);

    m=input('Do you want to change filters, y/n:','s');
    if m=='n'
        break
    end
end

%save workspace variable
outname = strcat(sampleName,'_filterV3_tot',num2str(cutTot),...
    '-LRrat',num2str(cutRat));
save(outname,'CTcut');
%dlmwrite(outfile 'near','Delimiter','\t');

%save and close figures
savefig(f1,figfile);
savefig(f2,outname);
close all;

%Generate OUTFILE that can be passed back to Goodman workflow
    %save as cell array 
%resize cell array, incorporate filtered data
chrID = chrID(1:aboveCut);
passcell = [chrID num2cell(CTcut(:,1:4))];
%write cell array to tab-delimted text
fileID = fopen(passfile,'w');
formatSpec = '%s\t %d\t %d\t %d\t %d\n';
[nrows,ncols] = size(passcell);
for row = 1:nrows
    fprintf(fileID,formatSpec,passcell{row,:});
end
fclose(fileID);

end


