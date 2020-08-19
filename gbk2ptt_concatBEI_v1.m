function gbk2ptt_concatBEI_v1(infile)
%last edit: 07-23-2020

%Takes .gbk of .gbff file and converts to .ptt (tab delimited)
%Concatenates contigs

%.ptt fields: 
%Location / Strand / Length / PID / Gene / Synonym / Code / COG / Product
%793..1353 / - / 561 / EIY83678.1 / - / 00001 / - / - / pooase protein

%Notes: 
    %"Synonym" field in .ptt is used for "GeneID" in 'mapped' file
    %'locus_tag' ("HMPREF1058_#####"), 'ID-'+5-digit number to --> 'Synonym'
    
%INPUT
    % Genbank annotation file with contigs
        %e.g. 'BEI.gbff'
%OUTPUT
    % .ptt table for concatenated BEI pseudo-genome
        %e.g. 'BEI.gbff.ptt'
%Note: ignore Warning for genbankread (line 333), irrelevant

    
%% Each contig put into an element of gbk structure
  %fields used for building output table
gbkS = genbankread(infile);
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
for i = 1:size(gbkS,2)
    s = gbkS(i);
    cds = s.CDS; %get CDSs for contig(i)
    len=size(cds,2);
    ci={cds.indices};
    %adjust indices for concatenation
    if i>1
    adj = (adj + str2double(gbkS(i-1).LocusSequenceLength));
    end
    %allocate cell arrays to store output
    loc = cell(len,1);
    siz=loc; orientation=loc; locusTag=loc; descr=loc;
    for j = 1:len
        a = ci{j}; %indices,[start stop]
        dif = a(2)-a(1);
        a = a + adj;
        start=num2str(a(1));
        stop=num2str(a(2));
        if dif > 0
            orientation(j) = {'+'};
            loc(j) = {strcat(start,'..',stop)};
        else
            orientation(j) = {'-'};
            loc(j) = {strcat(stop,'..',start)};
        end  
        siz(j) = {num2str(abs(dif))};
        %specific to BEI.gbff; changed 08-14-2020 to include 'ID-'
        locusTag(j) = {strcat('ID-',(cds(j).text(2,24:28)))};
        descr(j) = join(cellstr(cds(j).product)); 
    end
    %add Protein_ID & Gene
    pid = {cds.protein_id}';
    gene = {cds.gene}';
    %add '-' for Code & COG fields
    code = cell(len,1);
    code(:) = {'-'};
    cog = code;
   
%Location / Strand / Length / PID / Gene / Synonym / Code / COG / Product
%793..1353 / - / 561 / EIY83678.1 / - / 00001 / - / - / pooase protein
    %consolidate output in cell array
    itable = [loc, orientation, siz, pid, gene, locusTag,...
         code, cog, descr];
    table = [table; itable];
end

%% Output with column headers as tab delimited text
outfile = strcat(infile,'.ptt');
T = cell2table(table,'VariableNames',{'Location','Strand','Length',...
    'PID','Gene','Synonym','Code','COG','Product'});
writetable(T,outfile,'FileType','text','Delimiter','tab');
end


    