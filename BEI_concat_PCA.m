function BEI_concat_PCA(filterName,filterOn,norm_method,gbFile)
%last edit: August-29-2020
%requires natsortfiles.m (Matlab file exchange)
    
% Do PCA analysis on count-normalized data
    % normalization by:
        % CLR (centered log ratio) transformed data 
        % z-score
        
% BEI_concat_PCA('myTAsites',1,'clr','BEI_concat.gbff.txt')
% INPUTs
    % <filterName> is a string; sample name from file (eg. 'libraryTA')
    % <filterOn> is a flag: 1 for yes, 0 for...
        % 1== use Tn sites in filerName file, add pseudocount if absent
        % 0== find whole-sample-set intersect of TA-sites with non-zero entries
    % <norm_method> is a string: either 'clr' or 'zscore'  
    % <gbFile> is filename string; genbank formatted file for applying annotations to coeff
% OUTPUTs
    % plot of PC1 vs PC2 scores for all samples
    % plot of PC# vs %explained variation
    % annotated table of loadings
    % save workspace variables in .mat file
        

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
copyfile(filterFile, tempName);
clrtF = dlmread(tempName,'',0,1);
filter = sort(clrtF(:,1));
files=files(~idx); %remove filterFile from file list

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
ind = strfind(dirName,filesep);
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
  % Requires gbk file; limit to top-N insertion sites
  % Note: MUST maintain row index corresepondence of <coeff> to <TA site>

%Sort TA-site & coeff by absolute value of coeff
pc1_load = sortrows([z(1).data(:,1), coeff(:,1), abs(coeff(:,1))],3,'descend');
pc2_load = sortrows([z(1).data(:,1), coeff(:,2), abs(coeff(:,2))],3,'descend');
pc3_load = sortrows([z(1).data(:,1), coeff(:,3), abs(coeff(:,3))],3,'descend');
pc4_load = sortrows([z(1).data(:,1), coeff(:,4), abs(coeff(:,4))],3,'descend');

% Top N for PC1-PC4 in array
topN=50; %change as desired
pc_top_load = [pc1_load(1:topN,1:2); pc2_load(1:topN,1:2);...
    pc3_load(1:topN,1:2); pc4_load(1:topN,1:2)];

% Load genbank annotations
gbk = gbFile;
gbkS = genbankread(gbk);

%% Make row in table for each gbkS.CDS entry, loop through contigs
    %taken from gbk2table_concatBEI.m
table = {};
adj = 0;
for i = 1:length(gbkS)
    s = gbkS(i);
    cds = s.CDS; %get CDSs for contig(i)
    %collect contig info
    conNum = num2str(i);
    contig = cell(length(cds),1); contig(:) = {['>BEI_' conNum]};
    contigID = cell(length(cds),1); contigID(:) = {s.LocusName};
    if i>1
        adj = (adj + str2double(gbkS(i-1).LocusSequenceLength));
    end
    bounds = [cds.indices];
    bounds = [bounds(1:2:end);bounds(2:2:end)];
    bounds= bounds';%rows=CDS, columns=start,stop
    loc=sort(bounds,2,'ascend');%gene coordinates [low high]
    new_loc = loc + adj;
    first = num2cell(new_loc(:,1));
    last = num2cell(new_loc(:,2));
    size = num2cell(loc(:,2)-loc(:,1));
    prot = {cds.protein_id}';
    %grab the rest via loop
    orientation = cell(length(cds),1);
    geneID = cell(length(cds),1);
    descr = cell(length(cds),1);
    contig_indices = cell(length(cds),1);
    for j = 1:length(cds)
        if contains(cds(j).location,'complement')==1
            orientation(j) = {'-'};
        else
            orientation(j) = {'+'};
        end
        %specific to BEI.gbff; changed 08-14-2020 to include 'ID-'
        geneID(j) = {strcat('ID-',(cds(j).text(2,24:28)))};
        descr(j) = join(cellstr(cds(j).product));
        contig_indices(j) = {strcat(num2str(loc(j,1)),'..',num2str(loc(j,2)))};
    end
    %consolidate output in cell array
    itable = [contig, contigID, contig_indices...
        first, last, size, orientation, geneID, prot, descr];
    table = [table; itable];
end
genome_length = adj + str2double(s.LocusSequenceLength);

%% Annotate PC insertion site list using annotation table
coor = pc_top_load(:,1);
% Holders for: geneID, size, rel.position, product
len = length(coor);
%gene=cell(len,1); siz=[len,1]; pos=[len,1]; descr=cell(len,1);
%orientation=cell(len,1);
gene_info = cell(len,5);
loc=cell2mat(table(:,4:5));
for i = 1:len %for all coordinates in results table
    c = coor(i);%retrieve coordinate
    gidx = find(loc(:,1)<c,1,'last');%retrieve locus before insertion
   %intergenic handler
    if isempty(gidx)%insertion before first CDS
        name1 = cell2mat(table(1,8));
        geneID = {['pre_' name1]};
        siz = loc(1);%assumes first coordinate =1
        relPos = c/siz;
        descr = {'n/a'}; %table(1,10); %cell array
        orientation = {'+'};
    elseif c > loc(gidx,2) && gidx==length(table)%insertion after last CDS
        name1 = cell2mat(table(gidx,8));
        geneID = {['post_' name1]};
        siz = genome_length-loc(gidx,2);
        relPos = (c-loc(gidx,2))/siz;
        descr = {'n/a'}; %table(gidx,10);
        orientation = {'+'};
    elseif c > loc(gidx,2) && gidx~=length(table)%intergenic
        name1 = cell2mat(table(gidx,8));
        name2 = cell2mat(table(gidx+1,8));
        geneID = {[name1 '_' name2]};
        siz = loc(gidx+1,1)-loc(gidx,2);
        relPos = (c-loc(gidx,2))/siz;
        descr = {'n/a'}; %{cell2mat([table(gidx,10),{'_'},table(gidx+1,10)])};
        orientation = {'+'};
    else 
        %extract & store CDS info
        geneID = table(gidx,8);
        siz = loc(gidx,2)-loc(gidx,1);
        orientation = table(gidx,7);
        if cell2mat(orientation) == '+'
            relPos = (c-loc(gidx,1))/siz;
        else
            relPos = (loc(gidx,2)-c)/siz;
        end
        descr = table(gidx,10);
    end
    %consolidate output
    gene_info(i,:) = [geneID {siz} orientation {relPos} descr];
end

%% Export table
pc1=cell(topN,1); pc1(:)={'PC1'}; pc2=cell(topN,1); pc2(:)={'PC2'};
pc3=cell(topN,1); pc3(:)={'PC3'}; pc4=cell(topN,1); pc4(:)={'PC4'};
pcs=[pc1;pc2;pc3;pc4];

input = num2cell(pc_top_load); %coordinates and coeff
%csiz = num2cell(siz); cpos = num2cell(pos);
table_out = [pcs, input(:,2), input(:,1), gene_info];

Tout = cell2table(table_out,'VariableNames',{'PC','Coeff',...
    'TAsite','Gene','Size','Orientation','TnPosition','Product'});
outName = strcat('PCloading-',dirName,'_',sitesName,'_',norm_method);
writetable(Tout,outName,'FileType','text','Delimiter','tab');

end

%still ToDo:
    %use built in z-score ?
    %user-input for labeling sample groups ?
