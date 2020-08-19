%move genome name from end to sample, for contig ID

files = dir('INSEQ*processed*'); 
%reorder from UNIX to natural numeric sequence
files = natsortfiles({files.name});
%remove files in Goodman 'results' directory that aren't read tables
idx=~contains(files,'filter_cpm');
files=files(idx);

for i = 1:length(files)
    infile = char(files(i));
    %identify strings to modify
    sampleStart = strfind(infile,'scarf_') +6;
    sampleEnd = strfind(infile,'.bowtiemap') -1;
    sampleName = infile(sampleStart:sampleEnd);
    genomeStart = strfind(infile,'ed.txt_') +6;
    new_name = strcat(infile(1:sampleEnd),infile(genomeStart:end),...
        infile(sampleEnd+1:genomeStart-1));
    copyfile(infile,new_name);
end