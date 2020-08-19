function out = load_INseq_dataTable(sampleName)
%last edit: May-21-2020
%requires natsortfiles.m

%Imports numbers and saves them as workspace variable given by user
%INPUT
    % <sampleName> is just "SAMPLE" from Goodman file format:
      % "INSEQ_experiment.scarf_SAMPLE.bowtiemap_processed.txt_GENOME"
%OUT
    % first column is genome coordinate of Tn insertion (TA site)
    % second & third columns are #of reads mapped to L and R of insertion
    % fourth column is total #of reads mapped to insertion

%% action
%input list of file names, held in structure
files = dir('INSEQ*processed*'); 
%reorder from UNIX to natural numeric sequence
files = natsortfiles({files.name});
%remove files in Goodman 'results' directory that aren't read tables
files=files(~contains(files,'filter_cpm'));
idx=contains(files,strcat('.scarf_',sampleName,'.bowtiemap')); 
fileName=cell2mat(files(idx));
tempName = 'tempFile.txt';
copyfile(fileName, tempName);
out = dlmread(tempName,'',0,1);%column offset=1 to skip genomeName
out = sortrows(out,1);
delete(tempName);
end

%{
%for importing genome_TAsitesID
cell=readtable('Akk_TAsitesID.txt');
cell=table2cell(cell);
TAid=cell(:,2);
%}
