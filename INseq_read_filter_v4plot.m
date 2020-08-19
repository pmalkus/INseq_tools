function INseq_read_filter_v4plot(infile)
%last edit: April-24-2019

%Plot # unique insertions VS cutoff for total #reads per coordinate
    %L or R allowed to =0 (user input to change?)

% INFILE is read mapping table from from INseq pipeline: 
  %"INSEQ_experiment.scarf_<SAMPLE>.bowtiemap_processed.txt_<GENOME>"


%Initial/default cutoffs
cutTot = 1; %total number of reads (keep >= cutTot)
minLR = 0; %minimum number of L or R reads (keep >= minLR)

%Retrieve sample name for file output
sampleStart = strfind(infile,'scarf_') +6;
sampleEnd = strfind(infile,'.bowtie') -1;
sampleName = infile(sampleStart:sampleEnd);

%Temporarilty save original file with a new name understood by 'readtable'
ind1 = strfind(infile,'_experiment') +10;
tempName = strcat(infile(1:ind1),'.txt');
copyfile(infile, tempName);
%read numeric columns from infile, sort by insertion coordinate
clrt = dlmread(tempName,'','B1..E(end)');%change '' to '\t' for tab only
clrt = sortrows(clrt,1);
delete(tempName);

%%
function out = makeUnique_sum(in)
%Handle multiple entries for chr.coordinate (concatenated infile)
    %take sum of multiples
%and determine if multiples are present...
  %make_unique is very slow with large arrays, so skip if not needed
[u,ia,ib]=unique(in(:,1));
if length(u)==size(in,1)
    out = in;
else
    ul=[];ur=[];ut=[];
    for i=1:length(u)
        ul(i) = sum(in(ib==i,2));
        ur(i) = sum(in(ib==i,3));
        ut(i) = sum(in(ib==i,4));
    end
    ul=ul'; ur=ur'; ut=ut'; 
    uclrt = [u ul ur ut];
    uclrt = round(uclrt);
    out = uclrt;
end
clear u ia ib ul ur ut
end
%%
clrt = makeUnique_sum(clrt);

%Apply default filter and plot
CTcut = clrt((clrt(:,4) >= cutTot), :);
array=zeros(1000,2);
for j=1:1000
    above = clrt(clrt(:,4) >= j, :);
    array(j,:) = [j length(above)];
end
f1 = figure;
scatter(array(:,1),array(:,2));
set(gca,'Xscale','log');%,'Yscale','log'
set(gca, 'FontSize', 12)
name=strcat(sampleName,' (cutTotRead= ',num2str(cutTot),...
    ', minLR=',num2str(minLR),')');
title(name);
xlabel('cutoff for total # reads');
ylabel('#uniques insertions');
figfile = strcat(sampleName,'_filterV4_tot',num2str(cutTot),...
    '_minLR',num2str(minLR),'.fig');
%savefig(f1,figfile);

end


