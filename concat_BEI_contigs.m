function concat_BEI_contigs(sample_name,genome_name,contig_lengths)
%last edit: July-17-2020
%requires natsortfiles.m (from Matlab Central File Exchange)

% Concatenates read table contigs from sample files defined by user
% Call from current directory containing 'Goodman-format' output tables

%INPUT
    % <sample_name>, string for sample name
        % e.g. 'Pool1'
    % <genome_name>, string that comes before _contig# in filename
        % e.g. 'BEI'
    % <contig_lengths>, column vector of contig lengths, for contig #1-to-N)
        % e.g. [349378; 28485; 9345; 7798; 603]
%>> concat_BEI_contigs('Pool1','BEI',[349378; 28485; 9345; 7798; 603])

%OUTPUT
  % Concatented data table, same format as input tables
  % New sample name for outfile = genome_name'-concat' (e.g. BEI-concat)
  
%% Collect files and data
%input list of file names, held in structure
collect = strcat('INSEQ*processed.txt_',genome_name,'_*');
files = dir(collect); 
%reorder from UNIX to natural numeric sequence
files = natsortfiles({files.name});
%collect sample files
idx=contains(files,strcat('scarf_',sample_name,'.bowtiemap'));
files=files(idx);
%remove files that aren't read tables
idx=~contains(files,'filter_cpm');
files=files(idx);

%% Concatenate contigs
%Check that #of contigs matches # of sample files
len=length(contig_lengths);
if length(files)==len
    tempName='tempFile.txt';
    cat_data=[];
    adj=[0; contig_lengths(1:len-1)];
    for i=1:len
        infile = char(files(i));
        copyfile(infile, tempName);
        %read numeric columns from infile, sort by insertion coordinate
        clrt = dlmread(tempName,'','B1..E(end)');
        clrt = sortrows(clrt,1);
        %adjust coordinate numbering for contig
        clrt(:,1) = clrt(:,1) + sum(adj(1:i));
        cat_data = [cat_data; clrt];
    end
    delete(tempName);
else
    disp('SCRIPT TERMINATED!');
    disp('Number of contigs defined by user does not match #of sample files');
    return
end
%% Output consolidated data table in Goodman format
% Save output read table
passfile = strcat('INSEQ_experiment.scarf_',sample_name,...
    '.bowtiemap_processed.txt_',genome_name,'-concat.txt');
%resize cell array, incorporate filtered data
chrID={strcat('>',genome_name,'-concat')};
chrID(1:length(cat_data)) = chrID;
chrID = chrID';
passcell = [chrID num2cell(cat_data(:,1:4))];
%write cell array to tab-delimted text
fileID = fopen(passfile,'w');
formatSpec = '%s\t %d\t %d\t %d\t %d\n';
[nrows,ncols] = size(passcell);
for row = 1:nrows
    fprintf(fileID,formatSpec,passcell{row,:});
end
fclose(fileID);
end

