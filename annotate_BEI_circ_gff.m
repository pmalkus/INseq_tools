function annotate_BEI_circ_gff(gff_file,results_table)
%last edit: 08-12-2021

%Applies annotations to output from Goodman arrayed library well mapping
    % for fully assembled, circular Bv-BEI genome
    
%INPUTs
    % GFF annotation file (circular, one contig)
    % Arrayed match table
%OUTPUT
    % Annotations added to table and saved as new file
%PROCESS
    % For each Tn-insertion site from arrayed match table, find annotation
    % Intergenic insertions included as a separate class 

% Covert GFF to structure
gffObj = GFFAnnotation(gff_file);
gffStr = getData(gffObj);

% Read column 1,2 from results table and store as array
T = readtable(results_table, 'ReadVariableNames',false, 'Delimiter', 'tab');
map = table2array(T(2:end,1:2)); 
coor = str2double(map(:,2));
%to hold: geneID, class, size, insertion position (relative), product
gene_info = cell(length(coor),5);

% Find annotations
start = [gffStr.Start]';%not CDS start/stop, low-high coordinate bounds
stop = [gffStr.Stop]';
orientation = [gffStr.Strand]'; % can be: -/+/.
loc = double([start,stop]);%gene coordinates (location), [low high]
genL = 4878656; %Note: full genome length not in GFF file, taken from .gb

for j=1:length(coor)%for Tn-insertions in results_tale
    c = coor(j);%retrieve coordinate
    gidx = find(loc(:,1)<c,1,'last');%retrieve locus before insertion
    % Intergenic handler
    if isempty(gidx)% = insertion before first annotated gene
        ats = gffStr(1).Attributes;
        geneB = ats(strfind(ats,'.179.')+5:strfind(ats,';Name=')-1);
        geneID = {['pre_' geneB]};
        siz = loc(1,2);%assumes first coordinate =1
        relPos = c/siz;
        descr = 'n/a'; %ats(strfind(ats,';Name=')+6:end);
        class = 'intergenic';
    elseif c > loc(gidx,2) % = intergenic (ins.Coor larger than gene Stop)
        if gidx==length(loc)% = insertion after last annotated gene
            ats = gffStr(gidx).Attributes;
            geneA = ats(strfind(ats,'.179.')+5:strfind(ats,';Name=')-1);
            geneID = {['post_' geneA]};
            siz = genL-loc(gidx,2);
            relPos = (c-loc(gidx,2))/siz;
            descr = 'n/a'; %ats(strfind(ats,';Name=')+6:end);
            class = 'intergenic';
            gene_info(con_idx(j),:) = {geneID class siz relPos descr};
            continue
        end
  
        ats = gffStr(gidx).Attributes;
        geneA = ats(strfind(ats,'.179.')+5:strfind(ats,';Name=')-1);
        atsB = gffStr(gidx+1).Attributes;
        geneB = atsB(strfind(atsB,'.179.')+5:strfind(atsB,';Name=')-1);
        geneID = {['IG_' geneA '-' geneB]};
        siz = loc(gidx+1,1)-loc(gidx,2);
        relPos = (c-loc(gidx,2))/siz;
        descr = 'n/a';
        class = 'intergenic'; 
    else 
        %extract & store gene info
        ats = gffStr(gidx).Attributes;
        geneID = ats(strfind(ats,'.179.')+5:strfind(ats,';Name=')-1);
        siz = loc(gidx,2)-loc(gidx,1) +1;
        %realtive position of Tn-insertion requires orientation
        if orientation=='-'
            relPos = (loc(gidx,2)-c)/siz;
        else
            relPos = (c-loc(gidx,1) +1)/siz;
        end
        descr = ats(strfind(ats,';Name=')+6:end);
        class = gffStr(gidx).Feature;
    end
    %consolidate output
    gene_info((j),:) = {geneID class siz relPos descr};
end

%add to results table
tabl = table2array(T(2:end,1:end));
table_out = [tabl gene_info];

%output with column headers as tab delimited text
outfile = strcat(results_table(1:end-4),'_annotated');
Tout = cell2table(table_out,'VariableNames',{'Chromosome','Coordinate',...
    'StrainID','Plate','Well','Mismatches',...
    'GeneID','Class','Size','Position','Product'});
writetable(Tout,outfile,'FileType','text','Delimiter','tab');

end
    