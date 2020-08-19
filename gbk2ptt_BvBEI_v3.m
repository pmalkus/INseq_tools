function gbk2ptt_BvBEI_v3(infile)
%last edit, April-11-2019

%Takes .gbk of .gbff file and converts to .ptt (tab delimited)

%.ptt fields: 
%Location / Strand / Length / PID / Gene / Synonym / Code / COG / Product
%793..1353 / - / 561 / EIY83678.1 / - / 00001 / - / - / pooase protein

%Note: "Synonym" field in .ptt is used for "GeneID" in 'mapped' file

%BvBEI_v3: 
    %'locus_tag' ("HMPREF1058_#####"), 5-digit number only to --> 'Synonym'
    %'gene' to 'Gene' (dealing with lack of info)

%{
%for batch conversion from current directory
files = dir('*.gbff'); %input list of file names in structure
for i=1:length(files)
    infile = files(i).name;
    gbk2ptt_BvBEI_v3(infile);
end
%}
        
gbk = genbankread(infile);
cds=gbk.CDS;
c={cds.indices}';
len=max(size(cds));
ptt={};
%add range, strand, and length
for i=1:len
    a=c{i};
    dif = a(2)-a(1);
    start=num2str(a(1));
    stop=num2str(a(2));
    if dif > 0
        strand = '+';
        range = strcat(start,'..',stop);
    else strand = '-';
        range = strcat(stop,'..',start);
    end
    gene_length = num2str(abs(dif));
    ptt = [ptt; {range strand gene_length}];
end
%add Protein_ID
pid={cds.protein_id}';
ptt = [ptt pid];

%add gene (=Gene) locus_tag (=Synonym) and old_locus_tag (=Code)
    %for loop could break if no locus_tag found
        %make more resilient with independent searches
c={cds.text};
gene={}; synon={}; code={};
for i=1:len
    a=c{i};
    for j=1:size(a,1)
        if a(j,1:10)=='/locus_tag'
            %logical more robust (e.g 'contains')
            stin = strfind(a(j,:),'HMPREF1058_');
            syn = a(j,stin+11:stin+15); %locus_tag to Synonym field
        %else
            %syn = '-';
        end
    end
synon = [synon; {syn}];
end
%add '-' for COG field
fill = cell(len,1);
fill(:) = {'-'};
%'-' for Gene, Code, and COG fields; all empty in gbff for BvBEI
ptt = [ptt fill synon fill fill];

%add Product
 %most entries are on one line, but some are more - concatenate
c={cds.product}';
prod={};
for i=1:len
    out=c{i}(1,:);
    k=size(c{i},1);
    while j <= k
        mod = strcat({' '},c{i}(j,:));
        out = strcat(out,mod);
        j=j+1;
    end
    prod = [prod; out];
end
ptt = [ptt prod];

%output with column headers as tab delimited text
outfile = strcat(infile(1:strfind(infile,'.gb')),'ptt');
T = cell2table(ptt,'VariableNames',{'Location','Strand','Length',...
    'PID','Gene','Synonym','Code','COG','Product'});
writetable(T,outfile,'FileType','text','Delimiter','tab');


