function Tout = annotate_Akk_table_v3(fixGenome, coordinateFile, gbFile)
%last edit: July-14-2019
%works

% changes for v3:
  % added <fixGenome> flag for fixing known extra base in genome
    % sequence (Note: standard Goodman index file uses fixed genome sequence)
  % changed inputs to files; gbFile can be empty, in which case default is
  % used
  
%Applies annotations to input coordinates
%INPUT
    % <fixGenome>, flag for removing extra G at 1704819 in Akk genome
        % shifts downstream numbering from gbk for applying annotations
    % <coordinateFile>, TN insertion coordinate table (coordinates in 2nd column)
    % <gbFile>, annotation file (gbff or gbk)
%OUTPUT
    % Add annotations to table and save as new file
    % Handle intergenic insertions as a separate class
 %geneID / class / orientation / size / TN.position / product / notes
        %geneID = Amuc# (locus_tag); Amuc_####
        %class = CDS, tRNA, rRNA, ncRNA, pseudo, intergenic, misc
        %orientation = -, +, na
        %size = distance in bp btwn indices
        %TN.position = relative position of transposon
            %relative to orientation, when orientaton is given
            %relative to ascending coordinate when orientation 'na'
        %product = description of predicted product, na
        %notes = class specific information, e.g. sig-peptide for CDS
    
%PROCESS
    %For each CDS annotation query table and collect matching indices & data

%IMPROVE
    %Accommodate different types of tables: eg. -/+ header
    %Add "misc_binding" feature information = riboswitch (2x in Akk genome)
    
    
%%
%Format inputs  
   %NEEDs work to accommodate different types of tables: eg. -/+ header
infile=coordinateFile;
gbk=gbFile;
if isempty(gbk)
    gbk='/Users/pmalkus/Dropbox/Akk_etc_RVlab/INseq_things/Annotation_stuff/Amuc_CP001071.1.gb.txt';
    disp('using default Akkermansia genbank file, CP001071.1.gb.txt');
end
gbkStruct = genbankread(gbk);
f = featureparse(gbkStruct);

%Read column 2 from Goodman-style read table
tempName = 'tempFile.txt';
copyfile(infile, tempName);
coor=dlmread(tempName,'',0,1);%column offset=1 to skip genomeName
coor=coor(:,1);
coor=sort(coor);
delete(tempName);

%Holders for: geneID, class, orientation, size, rel.position,...
  %...product decription, signal peptide (-/+), notes
len = length(coor);
%annot = cell(length(coor),7);
gene=cell(len,1); class=cell(len,1); orient=cell(len,1);
siz=[len,1];pos=[len,1];
prod=cell(len,1);sigP=cell(len,1);note=cell(len,1);

%Collect all 'gene' info (assign pseudo); define intergenic regions
  %Use to cycle through other 'features' for matching 'product'(or 'note')
sq=f.gene;
indices=[sq.Indices];
%Fix indices from gbk file to match corrected genome (removed G at 1704819)
if fixGenome==1
    modi=find(indices < 1704819,1,'last');
    indices(modi+1:end)=indices(modi+1:end)-1;
end
indices=[indices(1:2:end);indices(2:2:end)];%gene coordinates, [start stop]
indices=indices';
loc=sort(indices,2,'ascend');%gene coordinates (location), [low high]  
for i=1:length(coor)%for all coordinates in results table
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
        if indices(gidx,2) > indices(gidx,1)
            orient(i) = {'+'};
        else
            orient(i) = {'-'};
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

%Signal peptide information
sig=f.sig_peptide;
siglocusID = {sig.locus_tag}';
for j=1:len
    g=gene(j);
    lidx=find(strcmp(siglocusID,g));
    if isempty(lidx)
        sigP(j) = {'-'};
    else
        sigP(j) = {'+'};
        note(j) = {sig(lidx).note};
    end
end
    
%Add annotations to results table
table = num2cell(coor);
csiz = num2cell(siz);
cpos = num2cell(pos);
out={}; out = [gene class orient csiz cpos prod sigP note];
table_out = [table out];

%NEED TO MAKE GENERAL, with respect to number of VariableNames (columns)
%output with column headers as tab delimited text
%outfile = strcat(results_table,'_annotated');
%Tout = cell2table(table_out,'VariableNames',{'Chromosome','Coordinate',...
%    'Plate','Gene','Class','Orientation','Size','Position','Product',...
%    'SignalP','Notes'}); %for array plate list
Tout = cell2table(table_out,'VariableNames',{'Coordinate',...
    'Gene','Class','Orientation','Size',...
    'Position','Product','SignalP','Notes'});%for AvB output table
%writetable(Tout,outfile,'FileType','text','Delimiter','tab');


end

%SCRAPS   
%[minDistance, indexOfMin] = min(abs(V-N));
%idx = find(strcmp([C{:}], 'bla')) %full string match
%idx = find(contains(C,'bla')) %partial string match accepted
    