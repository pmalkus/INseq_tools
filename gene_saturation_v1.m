function gene_saturation_v1(infile)

%Plot number of gene disrupted VS number unique insertions


%INPUT is Goodman mapped file
%INSEQ_experiment.scarf_SAMPLE.bowtiemap_processed.txt_GENOME_filter_cpm.txt_mapped"
    %Use aggregate of all data
    %find unique insertion coordinates
    %generate INPUT file of #insertion per gene
    
%Script does:
    %expand list of genes by number of insertions
    %sample randomly until all entries are done
    %collect #unique gene at each iteration
    %plot interation# vs #uniqueGene
    %(i.e. fraction of genes w/insertion -VS- number unique insertions sampled)
    
    
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
%randomize output
out = out(randperm(length(out)));

%build saturation array for plotting (gene identity not held)
newG=[1];
for i=2:length(out)
    if ismember(out(i),out(1:i-1)) == 1
        newG(i) = newG(i-1);
    else
        newG(i) = newG(i-1)+1;
    end
end

%save workspace variable
%Extract sample names & import data tables
sampleStart = strfind(infile,'scarf_') +6;
sampleEnd = strfind(infile,'.bowtie') -1;
sampleName = infile(sampleStart:sampleEnd);
outname = strcat(sampleName,'-geneSaturation-V1');
save(outname,'newG');

%Plotting
xcoord = [1:length(newC)]; %length(newC) should = length(out)
f1=figure;
plot(xcoord,newC);
set(gca, 'FontSize', 14);
title(outname,'FontSize', 16);
xlabel('number of insertion sites sampled','FontSize', 18);
ylabel('number of unique genes disrupted','FontSize', 18);
savefig(f1,outname);
end