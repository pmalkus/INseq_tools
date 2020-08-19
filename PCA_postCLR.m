function PCA_postCLR(filterName,filterOn,norm_method,gbFile)
%last edit: June-1-2020
%requires natsortfiles.m (Matlab file exchange)

%ToDo:
    %make annotation with genbank file optional ?
    %USE built in z-score ?
    %user-input for labeling sample groups ?
    
% Do PCA analysis on count-normalized data
    % normalization by:
        % CLR (centered log ratio) transformed data 
        % z-score
% user choices:
    % normalization method
    % handling zeros (TA sites present in Input but absent from Sample
        % user defined TA sites (filter) & pseduocount (when to apply?)
        % limit to coordinates in common for all inputs
        
% PCA_postCLR('BAPfeb19',1,'clr','Amuc_CP001071.1.gb.txt')
% INPUTs
    % <filterName> is a string; sample name from file (eg. 'Cecum1A')
    % <filterOn> is a flag: 1 for yes, 0 for...
        % 1== use Tn sites in filerName file, add pseudocount if absent
        % 0== find whole-sample-set intersect of TA-sites with non-zero entries
    % <norm_method> is a string: either 'clr' or 'zscore'  
    % <gbFile> is filename strring; genbank formatted file for applying annotations to coeff
% OUTPUTs
    % plot of PC1 vs PC2 scores for all samples
    % plot of PC# vs %explained variation
    % annotated table of loadings
        

%% Load all read tables in current directory
%input list of file names, held in structure
files = dir('INSEQ*processed*'); 
%reorder from UNIX to natural numeric sequence
files = natsortfiles({files.name});
%remove files in Goodman 'results' directory that aren't read tables
idx=~contains(files,'filter_cpm');
files=files(idx);

%% Filtering set-up, file handling
%filterString=[];
tempName = 'tempFile.txt';
idx=contains(files,filterName);
filterFile=char(files(idx));
files=files(~idx); %remove filterFile from file list
copyfile(filterFile, tempName);
clrtF = dlmread(tempName,'',0,1);
filter = sort(clrtF(:,1));

%if exist('filterName','var')
%Keep filter coordinates in vector outside stucture (as above)
%end

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
z=s; % use copy of structure for calculations, etc.

%% Define TA sites to use for analysis
    % Filter on user defined TA-sites, OR 
    % Reduce to TA sites present (non-zero counts) in all samples
lens=length(z);
if filterOn==1
    for i=1:lens
        %Filter based on filterFile
        z(i).data = z(i).data(ismember(z(i).data(:,1),filter),:);
        % Add back missing TA sites with 0
        inFilter_notData = setdiff(filter,z(i).data(:,1));
        fill = zeros(length(inFilter_notData),3);
        out = [z(i).data; inFilter_notData fill];
        z(i).data = sortrows(out,1);
    end
else % Reduce to TA sites present in all samples (all-samples intersect)
    %Find all-samples intersect of TA sites
    C=z(1).data(:,1);
    for i=1:lens-1
        [C, ~, ~]= intersect(C,z(i+1).data(:,1));
    end
    %Reduce data to all-samples intersect
    for i=1:lens
        %Reduce entries to sites in final entry (whole-set interesect)
        [~, ~, ib] = intersect(C,z(i).data(:,1));
        z(i).data = z(i).data(ib,:);
        z(i).data = sortrows(z(i).data,1); %necessary? instersect sorts?
    end
end

%% Normalize count data, user choice
norm_reads=[];
if strcmp(norm_method,'clr')==1
    % Extract reads, apply Centered Log Ratio transform
    for i=1:lens
        reads = z(i).data(:,4) +1; %add pseudocount to all before CLR
        z(i).clr = log2(reads/geomean(reads));
        norm_reads(:,i) = z(i).clr;
    end
%note: important that data in structure is sorted by coordinate!
    % rows of output matrix should be TA sites, columns are samples
elseif strcmp(norm_method,'zscore')==1
    % z-score: proportional normalization, log transform, scale to std-dev
    for i=1:lens 
        %Calculate scaling factor
        tot = sum(z(i).data(:,4));
        norm2cpm = 10^6/tot;
        %Apply scaling factor, add CPM data
        z(i).cpm = z(i).data(:,4) .* norm2cpm;
        cpm(:,i) = z(i).cpm;
    end
    cpm=cpm+1; %CPM pseudocount for all (added here so cpm=1 is log_cpm=0)
        %cpm_mean = mean(cpm,2);
        %cpm_std = std(cpm,0,2);
        %cpm_normVar = var(cpm,1,2);
    log_cpm = log(cpm); 
    log_cpm_mean = mean(log_cpm,2);
    log_cpm_std = std(log_cpm,0,2);
        %for std=0, replace with 1 to preserve 0 value
        %log_cpm_std(log_cpm_std==0) = 1;
    norm_reads=(log_cpm - repmat(log_cpm_mean,[1 lens])...
        ./ repmat(log_cpm_std,[1 lens]));
    % log_cpm_std=0 yields NaN for log_cpm=0 in all samples, replace with 0 
    norm_reads(isnan(norm_reads))=0;
end

% do PCA
[coeff,score,latent,~,explained,~]=pca(norm_reads');

% Save workspace variables; filename is name of current directory
dirName = pwd; 
ind = strfind(dirName,'/');
dirName = dirName(ind(end)+1:end);
save(dirName); 

%% Plotting
if filterOn==1
    sitesName=filterName;
else
    sitesName='nonZero-intersect';
end
% PCA scores & samples
figure; 
scatter(score(:,1),score(:,2),40, 'filled', 'r');
samples={s(:).sample}';
dx=1.5; dy=0;
text(score(:,1)+dx,score(:,2)+dy,samples,'FontSize',14);
%format, markers, labels, axes with explained values
xString = strcat('PC1 (',num2str(explained(1),3),'%)');
yString = strcat('PC2 (',num2str(explained(2),3),'%)');
xlabel(xString,'FontSize',16);
ylabel(yString,'FontSize',16);
title(strcat(dirName,': ',sitesName,',',norm_method),'FontSize',16);
%ax = gca;
%ax.YLabel.String = 'My y-Axis Label';
%ax.YLabel.FontSize = 14;
%ADD [user defined] color scheme for sample groups

% Explained -vs- PC
figure; bar(1:length(explained),explained);
xlabel('PC number','FontSize',16);
ylabel('% explained variation','FontSize',16);
title(strcat(dirName,': ',sitesName,',',norm_method),'FontSize',16);
xlim([0.5 10.5]); ylim([0 100]);
%ADD cumulative distribution 

%% Export annotated loadings (coeff)
  % Requires gbk file; limit to top N TA-sites, or top N%
  % Note: MUST maintain row index corresepondence of <coeff> to <TA site>
  %ADD 'Product' to output table ?!?!?

%Sort TA-site & coeff by absolute value of coeff
pc1_load = sortrows([z(1).data(:,1), coeff(:,1), abs(coeff(:,1))],3,'descend');
pc2_load = sortrows([z(1).data(:,1), coeff(:,2), abs(coeff(:,2))],3,'descend');
pc3_load = sortrows([z(1).data(:,1), coeff(:,3), abs(coeff(:,3))],3,'descend');
pc4_load = sortrows([z(1).data(:,1), coeff(:,4), abs(coeff(:,4))],3,'descend');

% Top N for PC1-PC4 in array
topN=50; %change as desired
pc_top_load = [pc1_load(1:topN,1:2); pc2_load(1:topN,1:2);...
    pc3_load(1:topN,1:2); pc4_load(1:topN,1:2)];

% Annotate (make optional based on user addition of gbkFile?)
coor = pc_top_load(:,1);
% Load genbank annotations
gbk = gbFile;
gbkStruct = genbankread(gbk);
f = featureparse(gbkStruct);
%Holders for: geneID, size, rel.position
len = length(coor);
gene=cell(len,1); siz=[len,1]; pos=[len,1]; prod=cell(len,1);
%Collect all 'gene' info (assign pseudo); define intergenic regions
sq=f.gene;
indices=[sq.Indices];%Note: genome numbering not adjusted for NCBI error
indices=[indices(1:2:end);indices(2:2:end)];%gene coordinates, [start stop]
indices=indices';
loc=sort(indices,2,'ascend');%gene coordinates (location), [low high]  
for i=1:len %for all coordinates in results table
    c=coor(i);%retrieve coordinate
    gidx=find(loc(:,1)<c,1,'last');%retrieve locus before insertion
    %Intergenic handler (assign all intergenic as orientation = '+')
    if isempty(gidx)%special case: insertion before first CDS
        name1 = {'beg--'}; name2 = {sq(1).locus_tag};
        name12 = [name1 name2];
        gene(i) = {strjoin(name12,'--')};
        siz(i)=loc(1);%assumes first coordinate =1
        pos(i) = c/siz(i);
        class(i) = {'intergenic'}; orient(i) = {'+'};
    elseif c > loc(gidx,2) %is integenic?
        if gidx==length(sq) %special case: insertion after last gene
            name1 = {sq(gidx).locus_tag}; name2 = {'--end'};
            name12 = [name1 name2];
            gene(i) = {strjoin(name12,'--')};
            siz(i)=f.source.Indices(2)-loc(gidx,2);
            pos(i) = (c-loc(gidx,2))/siz(i);
            class(i) = {'intergenic'}; orient(i) = {'+'};
            continue
        end
        name1 = {sq(gidx).locus_tag}; name2 = {sq(gidx+1).locus_tag};
        name12 = [name1 name2];
        gene(i) = {strjoin(name12,'--')};
        siz(i) = loc(gidx+1,1)-loc(gidx,2);
        pos(i) = (c-loc(gidx,2))/siz(i);
        class(i) = {'intergenic'}; orient(i) = {'+'}; 
    else 
        %for genes: extract & store info
        gene(i) = {sq(gidx).locus_tag};
        siz(i) = loc(gidx,2)-loc(gidx,1);
        pos(i) = (abs(c-indices(gidx,1)))/siz(i);
        if sq(gidx).pseudo==1
            class(i) = {'pseudo'};
        else
            class(i) = {'gene'};
        end
    end
end
siz=siz'; pos=pos';

%Make a table with relevant fields for sub-structure of Features (f)
  %class, locus_tag, product, note --> from f.tRNA, f.CDS, f.rRNA, f.ncRNA 
classID(1:length(f.tRNA),1) = {'tRNA'};
classID(end:end+length(f.CDS),1) = {'CDS'};
classID(end:end+length(f.rRNA),1) = {'rRNA'};
classID(end:end+length(f.ncRNA),1) = {'ncRNA'};
locusID = [{f.tRNA.locus_tag}';{f.CDS.locus_tag}';{f.rRNA.locus_tag}';{f.ncRNA.locus_tag}'];
prodID = [{f.tRNA.product}';{f.CDS.product}';{f.rRNA.product}';{f.ncRNA.note}'];
%Find annotation by locus_tag; add to holder cell arrays
    %product decription, signal peptide (-/+), notes
for j=1:len
    g=gene(j);
    lidx=find(strcmp(locusID,g));
    if isempty(lidx)
        prod(j)={'n/a'}; %NEED to deal with intergenics
        continue
    end
    prod(j)=prodID(lidx);
    if strcmp(class(j),{'gene'})==1
        class(j)=classID(lidx);
    end
end

%Export table
pc1=cell(topN,1); pc1(:)={'PC1'}; pc2=cell(topN,1); pc2(:)={'PC2'};
pc3=cell(topN,1); pc3(:)={'PC3'}; pc4=cell(topN,1); pc4(:)={'PC4'};
pcs=[pc1;pc2;pc3;pc4];

input = num2cell(pc_top_load); %coordinates and coeff
csiz = num2cell(siz);
cpos = num2cell(pos);
table_out = [pcs, input(:,2), input(:,1), gene, csiz, cpos, prod];

Tout = cell2table(table_out,'VariableNames',{'PC','Coeff',...
    'TAsite','Gene','Size','Position','Product'});
outName = strcat('PCloading-',dirName,'_',sitesName,'_',norm_method);
writetable(Tout,outName,'FileType','text','Delimiter','tab');

