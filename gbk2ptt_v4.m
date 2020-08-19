function gbk2ptt_v4(infile)

%take gbk file and convert to .ptt (tab delimited)
    %designed gbk from BAA-835 refseq
%.ptt fields: 
%Location / Strand / Length / PID / Gene / Synonym / Code / COG / Product
%793..1353 / - / 561 / WP_012419063.1 / yfgA / Amuc_0001 / AMUC_RS00010	/ - / pooase

%Version4: places RS##### where Amuc_#### is absent
    %NOTE - corrects coordinates for 3G-->2G error at 1704819
        %from output, manually fix errors in Amuc_1422 annotation
            %delete RS# (2x)  -->  RS12060, RS07610
            %insert Amuc_1422  -->  1704235..1707411, 3717 bp

gbk = genbankread(infile);
cds=gbk.CDS;
c={cds.indices}';
len=max(size(cds));
ptt={};

%add range, strand, and length
for i=1:len,
    a=c{i};
    if a(1) > 1704820  %modify indices to account for deleted base
        a = a-1;
    end
    dif = a(2)-a(1);
    low=num2str(a(1));
    high=num2str(a(2));
    if dif > 0,
        strand = '+';
        range = strcat(low,'..',high);
    else strand = '-';
        range = strcat(high,'..',low);
    end
    length = num2str(abs(dif));
    ptt = [ptt; {range strand length}];
end
%add Protein_ID
c={cds.protein_id}';
ptt = [ptt c];

%add gene (=Gene) old_locus_tag (=Synonym) and locus_tag (=Code)
    %for loop could break if no locus_tag found
        %make more resilient with independent searches
c={cds.text};
gene={}; synon={}; code={};
for i=1:len,
    a=c{i};
    for j=1:6,
        if a(j,1:6)=='/locus'
            cod = a(j,18:24); %locus_tag for Code field
            ji=j;
        end
    end
    if a(ji-1,1:5)=='/gene'
        ind = strfind(a(ji-1,:),'"');
        ind = [ind(1)+1:ind(2)-1];
        gen = a(ji-1,ind); %gene for Gene field
    else gen = '-';
    end
    if a(ji+1,1:4)=='/old'
        syn = a(ji+1,17:25); %old_locus_tag for Synonym field
    else syn = a(ji,18:24); %if absent, use locus_tag
    end
gene = [gene; {gen}];
synon = [synon; {syn}];
code = [code; {cod}];
end
ptt = [ptt gene synon code];

%add '-' for COG field
cog = cell(len,1);
cog(:) = {'-'};
ptt = [ptt cog];

%add Product
 %most entries are on one line, but some are more - concatenate
c={cds.product}';
prod={};
for i=1:len,
    out=c{i}(1,:);
    k=size(c{i},1);
     
    while j <= k,
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


