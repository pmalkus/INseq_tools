function gene_saturation_v1(infile)
%last edit: 12-01-2020

%Plots number of genes disrupted VS number unique insertions

%INPUT is Goodman mapped file
%"INSEQ_experiment.scarf_SAMPLE.bowtiemap_processed.txt_GENOME_filter_cpm.txt_mapped"
    
%PROCESS:
    %expand list of genes by number of insertions
    %sample randomly until all entries are done
    %collect #unique gene with insertion at each iteration
    %plot: number insertions sampled (iteration#) -VS- #genes w/insertion
    
%Make original file with a new name understood by 'readtable'
ind1 = strfind(infile,'mapped'); %just to be sure file is there
tempName = strcat(infile(1:16),'.txt');
copyfile(infile, tempName);
%read in data, convert to useable type
T = readtable(tempName, 'ReadVariableNames',false, 'Delimiter', 'tab');
delete(tempName);
data = table2array(T(2:end,1:2));
genes = data(:,1);
geneInd = [1:length(genes)];
geneInd = geneInd';
insNum = str2double(data(:,2));
gene_ins = [geneInd insNum];

%expand gene list by number of insertions; sample sequentially
out =[];
for i= 1:length(geneInd)
    if gene_ins(i,2)==0
        continue
    else 
        add = gene_ins(i,1) .* ones(gene_ins(i,2),1);
        out = [out; add];
    end
end

%Saturation analysis
  %loop and take mean to smooth effects of randomization step
%randomize output list of genes(*insertions)
newG = zeros(100,length(out));
for j=1:100
    rout = out(randperm(length(out)));%randomize order of <out>
    %build saturation array for plot (go through out list, add 1 for new gene)
    newG(j,1) = 1;
    for i=2:length(rout)
        if ismember(rout(i),rout(1:i-1)) == 1
            newG(j,i) = newG(j,i-1);
        else
            newG(j,i) = newG(j,i-1)+1;
        end
    end
mNewG = mean(newG,1); %take mean of 100 iterations of randomization
end

%Extract sample names & import data tables
sampleStart = strfind(infile,'scarf_') +6;
sampleEnd = strfind(infile,'.bowtie') -1;
sampleName = infile(sampleStart:sampleEnd);
outname = strcat(sampleName,'-geneSaturation-V1');
%save(outname,'newG'); %save workspace variable

%Plotting
xcoord = [1:length(mNewG)];
f1=figure;
plot(xcoord,mNewG);
set(gca, 'FontSize', 14);
title(outname,'FontSize', 16);
xlabel('number of insertion sites sampled','FontSize', 18);
ylabel('number of genes disrupted','FontSize', 18);
%savefig(f1,outname);
end