function annotate_BEI_circ_table(genbank_file,results_table)
%last edit: 08-09-2021

%Applies annotations to output from Goodman arrayed library well mapping
    % for fully assembled, circular BEI genome
%INPUT
    % Genbank annotation file (.gb or .gff) ???
    % Arrayed match table
%OUTPUT
    % Add annotations to table and save as new file
    % Handle intergenic insertions as a separate class 
%PROCESS
    % For each CDS annotation, query table and collect matching indices

%Format inputs
gbkStruct = genbankread(genbank_file);
%Read column 1,2 from results table and store as cell array
T = readtable(results_table, 'ReadVariableNames',false, 'Delimiter', 'tab');
map = table2array(T(2:end,1:2)); 
coor = str2double(map(:,2));

%to hold: geneID, size, insertion location (relative), product description(s)
gene_info = cell(length(coor),4);
%Find annotations
%{
for i=1:length(gbkStruct)
    conNum = num2str(i);
    contig={['>BEI_' conNum]};
    %retrieve indices in table matching contig(i)
    con_idx = find(ismember(map(:,1), contig));
%}    
% Use coordinates from results_tables to find gene from array
s=gbkStruct;
cds=s.CDS;%get CDSs for contig(i)
bounds=[cds.indices];
bounds=[bounds(1:2:end);bounds(2:2:end)];%gene coordinates, [start stop]
bounds=bounds';
loc=sort(bounds,2,'ascend');%gene coordinates (location), [low high]
    
for j=1:length(coor)%for all rows of [super]contig
    c=coor(j);%retrieve coordinate
    gidx=find(loc(:,1)<c,1,'last');%retrieve locus in before insertion
    %intergenic handler
    if isempty(gidx)%insertion before first CDS
        name1 = cds(1).text(2,24:28);
        geneID={['pre_' name1]};
        siz=loc(1);%assumes first coordinate =1
        relPos = c/siz;
        descr=cds(1).product(1,:);
    elseif c > loc(gidx,2) %intergenic (ins.Coord larger than CDS end)
        if gidx==length(cds)%insertion after last CDS
            name1 = cds(gidx).text(2,24:28);
            geneID={['post_' name1]};
            siz=str2num(s.LocusSequenceLength)-loc(gidx,2);
            relPos = (c-loc(gidx,2))/siz;
            descr=cds(gidx).product(1,:);
            gene_info(con_idx(j),:) = {geneID siz relPos descr};
            continue
        end
        name1 = cds(gidx).text(2,24:28);%BEWARE: fragile
        name2 = cds(gidx+1).text(2,23:28);
        geneID = {[name1 name2]};
        siz = loc(gidx+1,1)-loc(gidx,2);
        relPos = (c-loc(gidx,2))/siz;
        descr = {[cds(gidx).product(1,:),'_',cds(gidx+1).product(1,:)]};
    else 
        %extract & store CDS info
        geneID = cds(gidx).text(2,24:28);
        siz = loc(gidx,2)-loc(gidx,1);
        relPos = (abs(c-bounds(gidx,1)))/siz;
        descr = cds(gidx).product(1,:);
    end
    %consolidate output
    gene_info((j),:) = {geneID siz relPos descr};
end

%add to results table
tabl = table2array(T(2:end,1:end));
table_out = [tabl gene_info];

%output with column headers as tab delimited text
outfile = strcat(results_table,'_annotated');
Tout = cell2table(table_out,'VariableNames',{'Chromosome','Coordinate','StrainID',...
    'Plate','Well','Mismatches','GeneID','Size','Position','Description'});
writetable(Tout,outfile,'FileType','text','Delimiter','tab');

end

%SCRAPS   
%[minDistance, indexOfMin] = min(abs(V-N));
%idx = find(strcmp([C{:}], 'bla')) %full string match
%idx = find(contains(C,'bla')) %partial string match accepted
    