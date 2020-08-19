function INseq_read_filter_v4(infile)
%last edit, May-16-2020

%Derived from "INseq_read_filter_v3.m"
    %File naming conventions to match "pipeline_v3.pl"
    %Preserve Goodman naming and structure, with Matlab filter optional
    %Make new file to preserve original/unfiltered read table
 %_v4
    %Sum reads for replicate entries of a coordinate (concat. input file)
    %Plot # unique insertion found VS totRead cutoff
    %Apply cut-offs
    %L or R allowed to =0, user input to change
    
% Process primary data table output from INseq pipeline
    %Sum concatenated inputs
    %Apply cutoffs, plot, revise
    %Number of unique insertions, and matrix of coordinate/reads
    %Output file can be used for second-set of mapping-jobs (_SAMPLE.job2)

% INFILE is read mapping stats from from INseq pipeline: 
  %"INSEQ_experiment.scarf_<SAMPLE>.bowtiemap_processed.txt_<GENOME>"
  % preserve as "INSEQ_experiment.scarf_<SAMPLE>.bowtiemap_processed.txt_<GENOME>_RAW"
  % preserved original data table ("_RAW") can also be used as INFILE
    % Takes mean of concatenated input files
    % first column is genome name
    % second column is insertion position
    % third, fourth and fifth col = L, R, and Tot reads mapped to insertion
% OUTFILE is same format, saved as infile name:
  %"INSEQ_experiment.scarf_<SAMPLE>.bowtiemap_processed.txt_<GENOME>"

    
%Protection against overwriting RAW with filtered
inraw = strcat(infile,'_RAW');
files = ls;
if contains(files,inraw)==1
    disp('SCRIPT TERMINATED! Use RAW data file as input for filtering.');
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
minLR = 0; %minimum number of L or R reads (keep >= minLR)

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
% NB:only works with a single genome in sample
delete(tempName);

%%
function out = makeUnique_sum(in)
%Handle multiple entries for chr.coordinate (concatenated infile)
  %take sum of multiples, if present (speed improved, May 2020)   
[xnew,~,idx] = unique(in(:,1));
if length(xnew)==size(in,1)
    out = in;  
else
    Lnew = accumarray(idx(:),in(:,2));
    Rnew = accumarray(idx(:),in(:,3));
    Tnew = accumarray(idx(:),in(:,4)); 
    out = [xnew Lnew Rnew Tnew];
end
end
%%
clrt = makeUnique_sum(clrt);

%Apply default filter and plot
CTcut = clrt((clrt(:,4) >= cutTot), :);
array=zeros(1000,2);
for j=1:1000
    above = clrt(clrt(:,4) >= j, :);
    array(j,:) = [j length(above)];
end
f1 = figure;
scatter(array(:,1),array(:,2));
set(gca,'Xscale','log');%,'Yscale','log'
set(gca, 'FontSize', 12)
name=strcat(sampleName,' (cutTotRead= ',num2str(cutTot),...
    ', minLR=',num2str(minLR),')');
title(name);
xlabel('cutoff for total # reads)');
ylabel('#uniques insertions');
figfile = strcat(sampleName,'_filterV4_tot',num2str(cutTot),...
    '_minLR',num2str(minLR),'.fig');
savefig(f1,figfile);

%%
%Query user to revise cutTot and minLR values
  %Apply, then output coverage map, #of insertions passing filter
  %Save variable () to file for normalization, annotation, etc...
while(1)
    display('current value of cutTot='); cutTot,
    cutTot = input('Enter new cutoff for # of Total mapped reads per coordinate  ');
    display('current value of minLR='); minLR,
    minLR = input('Enter new cutoff for # of L or R mapped reads per coordinate  ');
%create new matrix ('CTcut') to preserve info in 'ccnzr'
    %CTcut = ccnzr((ccnzr(:,5) >= (1/cutRat)),:); 
    %CTcut = CTcut((CTcut(:,5) <= cutRat),:);  
    CTcut = clrt((clrt(:,4) >= cutTot), :);
    CTcut = CTcut((CTcut(:,2) >= minLR), :);
    CTcut = CTcut((CTcut(:,3) >= minLR), :);
    aboveCut=size(CTcut,1);
    msg=strcat('Number of insertion sites above cutoffs (totReads>',num2str(cutTot),...
        ', minLR=',num2str(minLR),')  =',num2str(aboveCut));
    display(msg);
%re-plot #unique vs cutTot for current cutoffs
    array=zeros(1000,2);
    for j=cutTot:1000
        above = CTcut(CTcut(:,4) >= j, :);
        array(j,:) = [j length(above)];
    end
    f2 = figure;
    scatter(array(:,1),array(:,2));
    set(gca,'Xscale','log');
    %xl=xlim; yl=ylim; %use later to preserve scale
    set(gca, 'FontSize', 12)
    name=strcat(sampleName,' (cutTotRead= ',num2str(cutTot),...
        ', minLR=',num2str(minLR),')');
    title(name);
    xlabel('cutoff for total # reads)');
    ylabel('#uniques insertions');
    figfile = strcat(sampleName,'_totRead',num2str(cutTot),...
        '_minLR',num2str(minLR),'.fig');
%User can choose to revise
    m=input('Do you want to change your mind, y/n:','s');
    if m=='n'
        break
    end
end

%save workspace variable
outname = strcat(sampleName,'_filterV4_tot',num2str(cutTot),...
    '_minLR',num2str(minLR));
save(outname,'CTcut');
%dlmwrite(outfile 'near','Delimiter','\t');

%close figures
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


