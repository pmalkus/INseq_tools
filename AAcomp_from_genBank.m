function AAcomp_from_genBank(infile,plotType)

%take genBank (.gbk) or fasta (.txt) file input, return global amino acid usage
%user inputs:
  %infile = string with file name to analyze (in current directory)
  %plotType = output plot type, accepts 'pie' OR 'bar' as input

len = length(infile);
type = infile((len-2):len);
if type == 'gbk'
    gbk = genbankread(infile);
    cds = gbk.CDS;
    c = {cds.translation}'; %cell array of all CDS
    %concatenate all amino acid seqeunces into a single character vector
    all_AA = [c{:}];
    %add genome name and NCBI accession number to plot
    accNum = gbk.LocusName;
    source = gbk.Source;
    name = strcat(source,' (',accNum,')'); 
elseif type == 'txt'
    cds = fastaread(infile);
    all_AA = cds.Sequence;
    %add fasta header to plot
    name = cds.Header;
end
msg = strcat(num2str(length(all_AA)),' amino acids analyzed');
display(msg);
%make a figure
fig = figure;
count = aacount(all_AA,'Chart',plotType);
title(name,'FontSize',16,'Interpreter', 'none')
%'interpreter' handles underscore in Accession#
%save output AA count in .mat file too?
end
