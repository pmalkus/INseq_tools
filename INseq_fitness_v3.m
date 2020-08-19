function INseq_fitness_v3(d,t1_sampleName,t2_sampleName,libraryTA,goodTA)
%last edit, July-30-2020

%Changes for version-3
    % fixed bug in consolidating data from annotation table
        % >> unique(locus,'stable') --> 'stable' to preserve TA-site order
    % calculate fitness AFTER consolidating counts into genes
        % v2 took mean of TA-fit to get gene-fit, low count calc. problem
    % include count data in output table to use for QC (for TA and gene)
        % for samples used in calculation of W, minimum and mean
    % save Matlab workspace
    
%INPUTs
  % d = number of bacteria at t2 / number of bacteria at t1
  % t#_name = SAMPLE string in sample file name
  % libraryTA = SAMPLE string in read table file used for TA-site filtering
  % goodTA = SAMPLE string in read table file used for TA-site filtering
    % "INSEQ_experiment.scarf_SAMPLE.bowtiemap_processed.txt_GENOME_filter_CPM.txt"
% >> INseq_fitness_v3(4, 'Mucin-t1', 'Mucin-t2', 'libTAsites', 'myTAsites','Akk.gb')
% NOTE: annotation file (*.gb*) must be in current directory

%OUTPUT
    % Cell array & Table of fitness values per TA site
    % Cell array & Table of fittness values per gene
    
%ToDo?
    %use count data for weighting locus-based mean
    %handling zeros... rates are >=bound. MARK to keep track
        %also in genes?
        %zeros in both t1 and t2, or specific to t2? or specific to either?
        
%% Load read tables in current directory
%input list of file names, held in structure
files = dir('INSEQ*processed*'); 
%reorder from UNIX to natural numeric sequence
files = natsortfiles({files.name});
%remove files in Goodman 'results' directory that aren't read tables
idx=~contains(files,'filter_cpm');
files=files(idx);

%% Filtering set-up, file handling
%changed so 'sample1' doesn't retrieve 'sample13', etc. (July-17-2020)
%filter for known TA-sites in library
tempName = 'tempFile.txt';
idx=contains(files,strcat('.scarf_',libraryTA,'.bowtiemap'));
filterFile=char(files(idx));
copyfile(filterFile, tempName);
clrtF = dlmread(tempName,'',0,1);
filter = sort(clrtF(:,1));
%secondary, post-calculation filter for good TA sites
idx=contains(files,strcat('.scarf_',goodTA,'.bowtiemap'));
filterFile=char(files(idx));
copyfile(filterFile, tempName);
clrtF = dlmread(tempName,'',0,1);
goodTA_filter = sort(clrtF(:,1));
%limit structure to user-defined sample files
idx1=contains(files,strcat('.scarf_',t1_sampleName,'.bowtiemap'));
idx2=contains(files,strcat('.scarf_',t2_sampleName,'.bowtiemap'));
idx=logical(idx1+idx2);
files=files(idx);

%% Consolidate sample names & read data into structure
s = struct('sample',{},'data',[]);
tempName = 'tempFile.txt';
for i=1:length(files)
    infile = char(files(i));
    copyfile(infile, tempName);
    s(i).sample = infile(strfind(infile,'scarf_')+6:strfind(infile,'.bowtiemap')-1);
    s(i).data = dlmread(tempName,'',0,1);%column offset=1 to skip genomeName
end
delete(tempName);

%% Define TA sites to use for analysis
    % Filter on user defined TA-sites
lens=length(s);
for i=1:lens
    %Filter based on filterFile
    s(i).data = s(i).data(ismember(s(i).data(:,1),filter),:);
    % Add back missing TA sites with 0
    [inFilter_notData, zidx] = setdiff(filter,s(i).data(:,1));
    fill = zeros(length(inFilter_notData),3);
    out = [s(i).data; inFilter_notData fill];
    s(i).data = sortrows(out,1);
    %keep indices for TA sites with counts=0 in sample data
    s(i).zidx = zidx;
end
counts=[];
%% Add pseusdocount & calculate TA site frequency
for i=1:lens 
    % Add pseudocount: 1 read for all sites
    counts(:,i) = s(i).data(:,4);
    counts(:,i) = counts(:,i) +1;
    %Calculate frequency
    tot = sum(counts(:,i));
    s(i).freq = counts(:,i) ./ tot;
end

%% Perform fitness calculation for each TA site
%Nt == frequency of clone in population = site counts / total counts
Nt1 = s(1).freq; 
Nt2 = s(2).freq; 
%For every coordinate calculate a fitness W, based on equation given by 
    %van Opijnen and Camilli, Nature 2009
W = zeros(length(Nt1),1);
for i=1:length(Nt1)
    numer = log(Nt2(i)*(d/Nt1(i)));
    denom = log((1-Nt2(i))*(d/(1-Nt1(i))));
    W(i) = numer/denom;
end

%% Output1 = W values per TA site, annotated
%Filter based on TA sites in 'goodTA' (goodTA could be same as libraryTA)
goodTA_fitness = W(ismember(filter,goodTA_filter));
%coordinates in goodTA_filter with pseudocount in sample (either t1 or t2)
zidx = unique([s(1).zidx; s(2).zidx]);
[~,goodTA_zidx,~] = intersect(goodTA_filter,filter(zidx));%zeros in goodTA
% Mean of counts from t1 and t2 to use for QC rank
counts = counts(ismember(filter,goodTA_filter),:);
mean_counts = mean(counts,2);
% Min of counts from t1 and t2 to use for QC rank (added July-30-2020)
min_counts = min(counts,[],2);
% Filter frequency data too
Nt1 = Nt1(ismember(filter,goodTA_filter));
Nt2 = Nt2(ismember(filter,goodTA_filter));

%% Annotate
% requires annotation file in current directory ('.gb')
table_out = annotate_Akk_TAsites(1,goodTA_filter); %flag fixes NCBI genome
taFit = [num2cell(min_counts), num2cell(mean_counts),...
    num2cell(goodTA_fitness), table_out];
%save output table
Tout = cell2table(taFit,'VariableNames',{'minCounts','meanCounts',...
    'Fitness','TAsite','Gene','Size','Position','Product'});
outName = strcat('Fitness.TA_d',num2str(d),'_',t1_sampleName,'_',...
    t2_sampleName,'_',goodTA,'.txt');
writetable(Tout,outName,'FileType','text','Delimiter','tab');

%% Output2 = W values per gene/locus
%retrieve count data for weighting mean calculation?
%number of TA with no counts per gene (zeros), or low counts??
locus = table_out(:,2);
[uloc,ia,~] = unique(locus,'stable');
gene_counts = zeros(length(ia),2);
for j=1:length(ia)-1
    range = (ia(j):ia(j+1)-1);
    numIns(j) = length(range);
    %sum frequency for all sites in locus
    gNt1(j) = sum(Nt1(range));
    gNt2(j) = sum(Nt2(range));
    %sum counts for all TA sites in locus, for QC
    gene_counts(j,:) = sum(counts(range,:),1);
end  
%special treatment for last entry
range = (ia(end):length(locus));
numIns(length(ia)) = length(range);  
gNt1(length(ia)) = sum(Nt1(range));
gNt2(length(ia)) = sum(Nt2(range));
gene_counts(length(ia),:) = sum(counts(range,:),1); 
%keep counts for condition for using in QC
mean_gene_counts = mean(gene_counts,2);
min_gene_counts = min(gene_counts,[],2); %added July-30-2020

%calculate fitness for each gene
gW = zeros(length(gNt1),1);
for i=1:length(gNt1)
    numer = log(gNt2(i)*(d/gNt1(i)));
    denom = log((1-gNt2(i))*(d/(1-gNt1(i))));
    gW(i) = numer/denom;
end

%Output table
geneFit = [num2cell(min_gene_counts), num2cell(mean_gene_counts),...
    num2cell(gW), num2cell(numIns)'];
geneFit = [geneFit, table_out(ia,2:end)]; %annotation

%save output table
Tout = cell2table(geneFit,'VariableNames',{'minCounts','meanCounts',...
    'Fitness','Inserts','Locus','Size','Position','Product'});
outName = strcat('Fitness.Locus_d',num2str(d),'_',t1_sampleName,'_',...
    t2_sampleName,'_',goodTA,'.txt');
writetable(Tout,outName,'FileType','text','Delimiter','tab');

outName = strcat('workspace_d',num2str(d),'_',t1_sampleName,'_',...
    t2_sampleName,'_',goodTA,'.mat');
save(outName);

