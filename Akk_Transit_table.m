function Akk_Transit_table(fix1422,gb_file)
%last edit: Mar-05-2021
%e.g. >> Akk_Transit_table(1,'CP001071.1.gb.txt')

% Generates annotation tables [for Transit]

% INPUTs
    % <fix1422>, TRUE/FALSE, remove extra G at 1704819 in BAA-835 genome
        % shift downstream genome numbering for applying annotations
    % <gb_file>, filename in current directory
        % works with 'genbank' file format (.gb, .gbk), NOT .gff
        % compatible with GenBank (CP001071) & RefSeq (NC010655) resource
        	% format detected based on input file header info
    
% OUTPUTs: tab delimited tables 
    % 'transitTable'
        % Fields required for Transit, plus addition column for Class
            % Use to sort, or to replace Product for non-coding regions
        % Length for CDS in aa (inc. stop codon), all other in nucleotides
     % 'fullTable'
        % Additional columns/fields provide more annotation information
        % Length in nucleotides for all entries
    % Note: LocusID for GenBank = 'Amuc_####' or 'Amuc_R####'
               % and for RefSeq = 'RS#####'
    % 'noIntergenic'
        % table of locus-pairs that lack intergenic region (overlap, abut)
   
%% Parse inputs
% Format inputs
gbkS = genbankread(gb_file);
f = featureparse(gbkS); %structure with CDS, tRNA, etc.
% evaluate annotation file format (GenBank vs. RefSeq)
if strcmp(gbkS.LocusName,'CP001071') == 1
    source = 'genbank';
elseif strcmp(gbkS.LocusName,'NC_010655') == 1
    source = 'refseq';
else
    disp('Unexpected annotation file, attempting with RefSeq format.');
    source = 'refseq';
end

% Use f.genes as starting point, then add other information from f.CLASS
%% ...for GenBank annotation file
if strcmp(source,'genbank') == 1
    % Populate fields
    indices = {f.gene(:).Indices}'; %parse for orientation and length later
    geneID = {f.gene(:).locus_tag}'; %cell array of all locus tags
    nloci = length(geneID);
    class = cell(nloci,1);
    product = cell(nloci,1);
    note = cell(nloci,1);
    sigPep = cell(nloci,1);
    % Add pseudogene info in 'product' and 'class'
    pseudo = {f.gene(:).pseudo}';
    pseudo = cellfun(@islogical,pseudo);
    fill = {'pseduogene'};
    class(pseudo) = fill; 
    product(pseudo) = fill;
    % Add CDS info (class, product, note)
    [~,ia,ib] = intersect(geneID,{f.CDS(:).locus_tag});
    class(ia) = {'CDS'};
    product(ia) = {f.CDS(ib).product}';
    note(ia) = {f.CDS(ib).note}';
    % Add signal peptide prediction info (genbank only)
    [~,ia,ib] = intersect(geneID,{f.sig_peptide(:).locus_tag}); 
    sigPep(ia) = {f.sig_peptide(ib).note}';
    % Add tRNA info (class, product)
    [~,ia,ib] = intersect(geneID,{f.tRNA(:).locus_tag});
    class(ia) = {'tRNA'};
    product(ia) = {f.tRNA(ib).product}';
    % Add rRNA info (class, product)
    [~,ia,ib] = intersect(geneID,{f.rRNA(:).locus_tag});
    class(ia) = {'rRNA'};
    product(ia) = {f.rRNA(ib).product}';
    % Add ncRNA info (class, product, note) - ONLY ONE, check indexing works!!!
    [~,ia,ib] = intersect(geneID,{f.ncRNA(:).locus_tag}); 
    class(ia) = {'ncRNA'};
    product(ia) = {f.ncRNA(ib).ncRNA_class}';
    note(ia) = {f.ncRNA(ib).note}';
end

%% ...for RefSeq annotation file
% Input data from RefSeq annotation file
if strcmp(source,'refseq') == 1
    % Populate fields
    indices = {f.gene(:).Indices}'; %parse for orientation and length later
    geneID = {f.gene(:).locus_tag}'; %cell array of all locus tags
    old_tag = {f.gene(:).old_locus_tag}';
    gene = {f.gene(:).gene}';
    nloci = length(geneID);
    class = cell(nloci,1);
    product = cell(nloci,1);
    % Add pseudogene info in 'product' and 'class'
    pseudo = {f.gene(:).pseudo}';
    pseudo = cellfun(@islogical,pseudo);
    fill = {'pseduogene'};
    class(pseudo) = fill; 
    product(pseudo) = fill;
    % Add CDS info (class, product, note)
    [~,ia,ib] = intersect(geneID,{f.CDS(:).locus_tag});
    class(ia) = {'CDS'};
    product(ia) = {f.CDS(ib).product}';
    % Add tRNA info (class, product)
    [~,ia,ib] = intersect(geneID,{f.tRNA(:).locus_tag});
    class(ia) = {'tRNA'};
    product(ia) = {f.tRNA(ib).product}';
    % Add rRNA info (class, product)
    [~,ia,ib] = intersect(geneID,{f.rRNA(:).locus_tag});
    class(ia) = {'rRNA'};
    product(ia) = {f.rRNA(ib).product}';
    % Add ncRNA info (class, product, note) - ONLY ONE, check indexing works!!!
    [~,ia,ib] = intersect(geneID,{f.ncRNA(:).locus_tag}); 
    class(ia) = {'ncRNA'};
    product(ia) = {f.ncRNA(ib).product}';
   % Add tmRNA info (class, product) (refseq only)
    [~,ia,ib] = intersect(geneID,{f.tmRNA(:).locus_tag});
    class(ia) = {'tmRNA'};
    product(ia) = {f.tmRNA(ib).product}'; 
end
   
%% Parse indices (optional: correct NCBI sequencing error)

% Using indices, generate Start, Stop, Orientation, Length (nucleotides)
start = zeros(nloci,1);
stop = zeros(nloci,1);
orientation = cell(nloci,1);
geneLength = zeros(nloci,1);

% Find orientation, reorder bounds for ascending start-stop
bounds = cell2mat(indices);
orilog = bounds(:,1) > bounds (:,2); % TRUE for - orientation
orientation(~orilog) = {'+'};
orientation(orilog) = {'-'};
bounds(orilog,:) = [bounds(orilog,2),bounds(orilog,1)];% all [low, high]
bitend = f.source.Indices(2); %genome length

% Adjust numbering for NCBI sequencing error, fix annotation for Amuc_1422
if fix1422 == 1  %flag, from user input
    adj = find(bounds(:,1) > 1704819);
    bounds(adj,1) = bounds(adj,1)-1;
    adj = find(bounds(:,2) > 1704819);
    bounds(adj,2) = bounds(adj,2)-1;
    bitend = bitend-1; % adjust full genome length
    % edit entry for 1422
    if strcmp(source,'genbank') == 1
        ind = find(strcmp(geneID,'Amuc_1422'));
    elseif strcmp(source,'refseq') == 1
        ind = find(strcmp(old_tag,'Amuc_1422'));
    end
    % Define indices, etc.
    bounds(ind,:) = [1704235 1707411];
    orientation(ind) = {'-'};
    class(ind) = {'CDS'};
    product(ind) = {'recG'};
    note(ind) ={'added back after correcting NCBI sequencing error'};
end
geneLength = bounds(:,2) - bounds(:,1) +1;

% Convert to cell array
start = num2cell(bounds(:,1));
stop = num2cell(bounds(:,2));
geneLength = num2cell(geneLength);

% Check for gene-pairs without intergenic region
ino = []; %holder of indices(i) in geneID that have overlapping gene(i+1)
for i=1:nloci-1
    if bounds(i,2)+1 > bounds(i+1,1)-1
        %disp(['Gene-pairs without intergenic region: ',...
            %cell2mat(geneID(i)),' and ',cell2mat(geneID(i+1))]);
        ino = [ino i];
    end
end
% Save list of gene-pairs w/out intergenic for export
noInter = [geneID(ino) geneID(ino+1)];

%% Add intergenics
% Intergenics need to be added with the following fields
    % start / stop / orientation / gene_length / geneID / class / product
    
% Collect intergenics
%interL = ~ismember(1:nloci,ino); %logical for intergenics
interi = setdiff(1:nloci,ino);
interN = length(geneID) - length(ino); %number of intergenics
% Special case for beginning of genome
if bounds(1,1) >2
    istart(1) = 1;
    istop(1) = bounds(1,1)-1;
    igeneID(1) = strcat({'Beg--'},geneID(1));
end
% Continue with intergenics (none for overlapping genes)
istart(2:interN) = bounds(interi(1:end-1),2) +1;
istop(2:interN) = bounds(interi(1:end-1)+1,1)-1;
igeneID(2:interN) = strcat(geneID(interi(1:end-1)),...
    {'--'},geneID(interi(1:end-1)+1));
% Special case for end of genome
if bounds(end,2) < bitend-2
    istart(length(istart)+1) = bounds(end,2)+1;
    istop(length(istop)+1) = bitend;
    igeneID(length(igeneID)+1) = strcat(geneID(end),{'--End'});
end
ilength = istop - istart +1;

% Add orientation, class, product (same for all intergenic)
iorientation(1:length(istart),1) = {'+'}; % all ascending
iclass(1:length(istart),1) = {'intergenic'};
iproduct(1:length(istart),1) = {'intergenic'};
% Invert vectors, convert to cell
istart = num2cell(istart)';
istop = num2cell(istop)';
ilength = num2cell(ilength)';
igeneID = igeneID';

% Append intergenic to table
start = [start; istart];
stop = [stop; istop];
orientation = [orientation; iorientation];
geneLength = [geneLength; ilength];
geneID = [geneID; igeneID];
class = [class; iclass];
product = [product; iproduct];
% Add length to unused fields
extra = nloci:length(start);
note(extra) = {'-'};
sigPep(extra) = {'-'};
old_tag(extra) = {'-'};
gene(extra) = {'-'};

%% Consolidate and export table
% For Transit output, calculate length in aa for CDS
cds_log = strcmp(class,'CDS');
cds_length = geneLength(cds_log);
cds_length = cell2mat(cds_length);
cds_aa = cds_length ./ 3;
cds_aa = num2cell(cds_aa);
transit_length = geneLength;
transit_length(cds_log) = cds_aa;

% Note adjustment for sequncing error in Amuc_1422 in filename of output
fixName = '';
if fix1422 ==1
    fixName = '_fixed1422';
end

% Output table for Transit (+Class column, delete or use to replace Product) 
empty(1:length(geneID),1)={'-'}; %filler column for transit
table_out = [product, start, stop, orientation, transit_length,...
    empty, empty, empty, geneID, class];
% Fill in all empty fields with '-'
log_array = cellfun(@isempty,table_out);
table_out(log_array) = {'-'};
% Output with column headers as tab delimited text
outfile = strcat(gb_file,fixName,'_transitTable.txt');
Tout = cell2table(table_out,'VariableNames',{'Product','Start','Stop',...
    'Orientation','Length','Empty1','Empty2','Empty3','LocusID','Class'});
writetable(Tout,outfile,'FileType','text','Delimiter','tab');

% Output full table from GenBank
if strcmp(source,'genbank') == 1
    table_out = [start stop orientation geneLength geneID class product note sigPep];
    % Fill in all empty fields with '-'
    log_array = cellfun(@isempty,table_out);
    table_out(log_array) = {'-'};
    % Output with column headers as tab delimited text
    outfile = strcat(gb_file,fixName,'_fullTable.txt');
    Tout = cell2table(table_out,'VariableNames',{'Start','Stop','Orientation',...
        'Length','LocusID','Class','Product','Note','SignalPeptide'});
    writetable(Tout,outfile,'FileType','text','Delimiter','tab');
end
% Output full table from RefSeq
if strcmp(source,'refseq') == 1
     table_out = [start stop orientation geneLength geneID class product old_tag gene];
    % Fill in all empty fields with '-'
    log_array = cellfun(@isempty,table_out);
    table_out(log_array) = {'-'};
    % Output with column headers as tab delimited text
    outfile = strcat(gb_file,fixName,'_fullTable.txt');
    Tout = cell2table(table_out,'VariableNames',{'Start','Stop','Orientation',...
        'Length','LocusID','Class','Product','OldLocusID','Gene'});
    writetable(Tout,outfile,'FileType','text','Delimiter','tab');
end   

% Output list of gene-pairs w/out intergenic region
outfile = strcat(gb_file,fixName,'_noIntergenic.txt');
Tout = cell2table(noInter,'VariableNames',{'Locus1','Locus2'});
writetable(Tout,outfile,'FileType','text','Delimiter','tab');

end