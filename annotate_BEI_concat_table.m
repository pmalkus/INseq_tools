function annotate_BEI_concat_table(gbk_file,results_table)
%last edit: 08-21-2020

%Applies annotations to output from Goodman arrayed library well mapping
    %for concatenated BEI genome
    
%INPUT
    % Annotation file in gbk or gbff format
    % Arrayed_match table
%OUTPUT
    % Add annotations to arrayed_match table and save as new file
        % Intergenic insertions handled as a separate class 
        % Note: Contig info not included in output, but could be
%PROCESS
    %Each contig put into an element of gbk structure
    %Build concatenated genomeo annotation table
    %For each insertion site coordinate in Array_match_table,
        %query table and collect matching annotation
%Note: ignore Warning for genbankread (line 333), irrelevant

%% Build annotation table for concatenated genome
%Format inputs
gbkS = genbankread(gbk_file);
%Read column 1,2 from results table and store as cell array
T = readtable(results_table, 'ReadVariableNames',false, 'Delimiter', 'tab');
map = table2array(T(2:end,1:2)); 
coor = str2double(map(:,2));

% Make row in table for each gbkS.CDS entry, loop through contigs
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
len = length(coor);
gene_info = cell(len,5);
loc = cell2mat(table(:,4:5));
for i = 1:len %for all coordinates in results table
    c = coor(i);%retrieve coordinate
    gidx = find(loc(:,1)<c,1,'last');%retrieve locus before insertion
   %intergenic handler
    if isempty(gidx)%insertion before first CDS
        name1 = cell2mat(table(1,8));
        geneID = {['pre_' name1]};
        siz = loc(1);%assumes first coordinate =1
        relPos = c/siz;
        descr = {'n/a'}; %table(1,10);
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

%% Add annotations to results table
table_out = table2array(T(2:end,1:end));
table_out = [table_out gene_info];

%output with column headers as tab delimited text
outfile = strcat(results_table,'_annotated');
Tout = cell2table(table_out,'VariableNames',{'Chromosome','Coordinate','StrainID',...
    'Plate','Well','Mismatches','GeneID','Size','Orientation','TnPosition','Product'});
writetable(Tout,outfile,'FileType','text','Delimiter','tab');

end

    