function gbk2table_concatBEI_v1(genbank_file)
%last edit: 07-22-2020

%Generates full annotation table from gbk/gbff with contigs
%INPUT
    % Genbank annotation file with contigs
        %e.g. 'BEI.gbff'
%OUTPUT
    % Annotation table for whole genome
    % Includes contig and concatenated genome indices
        %e.g. 'BEI.gbff_table.txt'
%Note: Contitg number assigned by order listed in gbk file (descending size)
%Note: ignore Warning for genbankread (line 333), irrelevant

    
%% Each contig put into an element of gbk structure
  %fields used for building output table
gbkS = genbankread(genbank_file);
%gbkS.LocusName == contig ID
%gbkS.LocusSequenceLength == contig length
%gbkS.CDS == geneInfo
        %.indices
        %.protein_id
        %.product
        %.text with /subfields
            %e.g. /locus_tag="HMPREF1058_03898"

%% Make row in table for each gbkS.CDS entry, loop through contigs
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

%% Output with column headers as tab delimited text
outfile = strcat(genbank_file,'_table.txt');
Tout = cell2table(table,'VariableNames',...
    {'Contig','ContigID','ContigIndices'...
    'First','Last','Size','Orientation',...
    'LocusTag','ProteinID','Description'});
writetable(Tout,outfile,'FileType','text','Delimiter','tab');

end
    