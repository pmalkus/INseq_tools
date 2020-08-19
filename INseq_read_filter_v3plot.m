function INseq_read_filter_v3plot(infile)
%Last edit: June-4-2019 (just changed cutRat default to 10^6)

%Plot L/R ratio VS total #reads per coordinate
    %Note: coordinates with either L or R =0 are filtered out

% INFILE is read mapping stats from from INseq pipeline: 
  %"INSEQ_experiment.scarf_<SAMPLE>.bowtiemap_processed.txt_<GENOME>"

%Initial/default cutoffs
cutTot = 1; %total number of reads (keep >= cutTot)
cutRat = 10^7; %ratio of L/R reads (keep <= cutRat)
minLR = 1; %number of L or R reads (keep >= minLR)
           %has to be >0 for ratio calculation!

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

%FILTERING
%remove rows with totRead < initial cutoff for total reads (=1)
carCut = clrt(clrt(:,4) >= cutTot, :);
%remove rows with L or R < default cutoff (+1)
carCut_noz = carCut((carCut(:,2)>=minLR),:);
carCut_noz = carCut_noz((carCut_noz(:,3)>=minLR),:);
%calculate ratio of L-to-R, apply default cutoff
ratLR=carCut_noz(:,2)./carCut_noz(:,3);
ccnzr = [carCut_noz ratLR];
ccnzr = ccnzr((ccnzr(:,5)>=(1/cutRat)),:);
ccnzr = ccnzr((ccnzr(:,5)<=(cutRat)),:);

%Plot total reads vs. ratio of LtoR, using DEFAULT CUTOFFs
f1 = figure;
scatter(ccnzr(:,4),ccnzr(:,5));
set(gca,'Xscale','log','Yscale','log');
xl=xlim; yl=ylim; %use later to preserve scale
set(gca, 'FontSize',14)
name=strcat(sampleName,' (cutTotRead= ',num2str(cutTot),...
    ', cutLRratio=',num2str(cutRat),')');
title(name, 'FontSize',20);
xlabel('total mapped reads (per coordinate)', 'FontSize',16);
ylabel('ratio of L-to-R mapped reads (per coordinate)', 'FontSize',16);
figfile = strcat(sampleName,'_filterV3_tot',num2str(cutTot),...
    '_LRrat',num2str(cutRat),'.fig');

end


