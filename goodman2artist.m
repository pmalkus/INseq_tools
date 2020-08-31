function [Control,Experiment] = goodman2artist(Input,Sample,TAsites)
%last edit: August-29-2020
%>> [Control,Experiment] = goodman2artist('Input1','69Cec',TAsites);

%Converts Goodman read table into Artist-ready workspace variables
  % Input data is filtered by TA sites; and TA sites found in the 
  % genome sequence, but not the data, are added with value=0.
%INPUTs
  % 'Input' is sample string for Goodman read table Input
  % 'Sample' is sample string for Goodman read table Test Sample
  %  <TAsites> variable must be in workspace (vector of TA sites in genome)
%OUTPUTS
  % <Control> = Total reads for Input sample
  % <Experiment> = Total reads for Test sample

%Create <TAsites> from Genome_TAsites.txt file
%TAsites_file=strcat(genome,'_TAsites.txt');
%T=readtable(TAsites_file);
%TAsites=table2array(T(:,2));

%Retrieve files
files = dir('INSEQ*processed*'); %input list of file names in structure
%Reorder from UNIX numeric sequence to standard
files = natsortfiles({files.name});
%remove files in Goodman 'results' directory that aren't read tables
idx=~contains(files,'filter_cpm');
files=files(idx);

%Use TAsites to fill & filter
[Control] = process_reads(Input);
[Experiment] = process_reads(Sample);

%% extract, filter and fill data
function [out] = process_reads(name)
%Load data tables from input file
idy=contains(files,strcat('.scarf_',name,'.bowtiemap')); 
infile=char(files(idy));
%Extract data tables
tempName = 'tempFile.txt';
copyfile(infile, tempName);
data = dlmread(tempName,'',0,1);%column offset=1 to skip genomeName
data = sortrows(data,1);
data = data(:,[1 4]);%keep only coordinate and total reads

%Remove non-TA sites from data file
data=data(ismember(data(:,1),TAsites),:);

%If absent from TAsites, add coordinate to data with read count =0
inTAnotData = setdiff(TAsites,data(:,1));
fill = zeros(length(inTAnotData),1);
data = [data; inTAnotData fill];
out = sortrows(data,1);
out = out(:,2);%reads for all TA sites

delete(tempName);
end
%% finish
end