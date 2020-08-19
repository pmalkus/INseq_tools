function Artist_run_directory(genome,inputName,filterName)
%last edit: May-13-2020, works

% Acts on all Test Samples in currennt directory

% Consolidates ARTIST workflow
% Generate normalized Input (mean of ControlSims)

% Parse all inputs and generate .mat file containing all workspace
    % variables needed for conditional gene analysis with ARTIST 
    
% INPUTs
    % 'Genome', prefix for Genome_TAsites and Genome_TAsitesID .txt files
    % 'inputName', filname: sample name for read table from Goodman pipeline
        % coordinates in column-2, total reads in column-5
    % 'sampleName', filname: sample name for read table from Goodman pipeline
        % coordinates in column-2, total reads in column-5
    % OPTIONAL: 'filterName', sample name for read table from Goodman pipeline
        % used to filter sample data prior to filling missing TA sites
        % leave out of command to use input data table as is
% OUTPUT
    % .mat file named by inputName_sampleName_filteName
    % contains workspace variables for analysis and MWU stats
        % also TAhits: Tn insertions per locus
        % and NormInput: mean of simulated/normalized counts for each site

        
% Loop through all Sample files in directory, exclude Input and Filter
    % Carefully define what Input and Filter shouold be
    % Include a single file for each in the current directory

    
%input list of file names, held in structure
files = dir('INSEQ*processed*'); 
%reorder from UNIX to natural numeric sequence
files = natsortfiles({files.name});
%remove Input and Filter files
idx=~contains(files,inputName);
files=files(idx);
idx=~contains(files,filterName);
files=files(idx);
%loop through all sample files
for i=1:length(files)
    infile = char(files(i));
    sampleName = infile(strfind(infile,'scarf_')+6:strfind(infile,'.bowtiemap')-1);
    Artist_run_v2(genome,inputName,sampleName,filterName);
end


