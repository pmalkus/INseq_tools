function table_out = annotate_Akk_TAsites(fixGenome,TAsites)
%last edit: July-30-2020

% Annotates column vector of genome coordinates, TA sites
% Output is cell array of coordinates and annotations
% Genbank file found in current dirctory from '.gb' (not user input)

%INPUT
    % <fixGenome>, flag for removing extra G at 1704819 in Akk genome
        % shifts downstream numbering from gbk for applying annotations
    % <TAsites> column vector of genome coordinates
        
coor = TAsites;
% Load genbank annotations
gbk = dir('*.gb*'); %find gbk file in curdir
gbk = gbk.name;
gbkStruct = genbankread(gbk);
f = featureparse(gbkStruct);
%Holders for: geneID, size, rel.position
len = length(coor);
gene=cell(len,1); siz=[len,1]; pos=[len,1]; prod=cell(len,1);
%Collect all 'gene' info (assign pseudo); define intergenic regions
sq=f.gene;
indices=[sq.Indices];
%For annotation, option to ADJUST indices for NCBI error in Amuc_1422
if fixGenome==1    
    modi=find(indices < 1704819,1,'last');
    indices(modi+1:end)=indices(modi+1:end)-1;
end 
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
%Output cell array for building table
ccoor = num2cell(coor);
csiz = num2cell(siz);
cpos = num2cell(pos);
table_out = [ccoor, gene, csiz, cpos, prod]; 
end