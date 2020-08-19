function INseq_read_filterFile_v2(lowFilter,infile,filterFile)
%last edit, May-16-2020

%version 2: 
    % add 'lowFilter' user input, flag for low number filtering (1==yes)
    % change sequence of filtering steps
%Overview
    %Filters "infile" data for noise: keep coordinates with totReads >= 3
    %Filters "infile" data table using insertion coordinates of "filterFile"
    %File naming conventions to match "pipeline_v3.pl"
    %Make new file to preserve original/unfiltered read table ("_RAW")
    %SUM reads for replicate entries of coordinates (concatenated input file)
    %Output file can be used for second-set of mapping-jobs (_SAMPLE.job2)
% OUTFILE is same format, and saved with infile name!!!
  %"INSEQ_experiment.scarf_<SAMPLE>.bowtiemap_processed.txt_<GENOME>"

%Changes to come?
    %make total read filter # a user input?
    %user input to choose sum or mean?
%%    
%Protection against overwriting RAW with filtered
inraw = strcat(infile,'_RAW');
files = ls;
if contains(files,inraw)==1
    display('SCRIPT TERMINATED! Use RAW data file as input for filtering.');
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

%
%% Retrieve data from infiles
function [sampleName, out_data, chrID] = import(input)
%Extract sample names & import data tables
sampleStart = strfind(input,'scarf_') +6;
sampleEnd = strfind(input,'.bowtie') -1;
sampleName = input(sampleStart:sampleEnd);
%Replace underscore in sampleName with dash
und = strfind(sampleName,'_');
if ~isempty(und)
    sampleName(und)='-';
end
%Collect input data in matrices
%Temporarilty save original file with a new name understood by 'readtable'
tempName = strcat(input(1:5),'.txt');
copyfile(input, tempName);
out_data = dlmread(tempName,'','B1..E(end)');
out_data = sortrows(out_data,1); %necessary?

%read character column from infile and store as cell array for output
T = readtable(tempName, 'ReadVariableNames',false, 'Delimiter', 'tab');
chrID = table2array(T(:,1)); %resize before output
% NB:only works with a single genome in sample

delete(tempName)
end
%%
[sample, clrt, chrID] = import(infile);
[~, clrtF, ~] = import(filterFile);

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

if lowFilter==1
%Apply noise filter (before summing replicates!)
totReads = 3; %NOTE: Goodman script uses >3 total reads filter (not >=3)
clrt = clrt((clrt(:,4) >= totReads), :);
end

%Apply filterFile: keep coordinates of infile that are in filterFile
filter = clrtF(:,1);
tf = ismember(clrt(:,1),filter);
fclrt = clrt(tf,:);

%Combine multiple entries for a single coordinate (sum)
    %do after other filtering steps for speed
fuclrt = makeUnique_sum(fclrt);
out = sortrows(fuclrt,1);

%Generate OUTFILE that can be passed back to Goodman workflow
    %save as cell array 
%resize cell array, incorporate filtered data
chrID = chrID(1:size(out,1));
passcell = [chrID num2cell(out)];
%write cell array to tab-delimted text
fileID = fopen(passfile,'w');
formatSpec = '%s\t %d\t %d\t %d\t %d\n';
[nrows,ncols] = size(passcell);
for row = 1:nrows
    fprintf(fileID,formatSpec,passcell{row,:});
end
fclose(fileID);
end


