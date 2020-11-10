function INseq_AvsB_filterF_v3_plot(mark, sampleA, sampleB, filterF)
% last edit, Nov-07-2020
%>> INseq_AvsB_filterF_v3_plot([],'GF-1','GF-2','lib-Tn')
%>> INseq_AvsB_filterF_v3_plot([123, 1234, 12345],'GF-1','GF-2','lib-Tn')
  
% v3 changes from v2:
    % updated file name input (just sample name, not whole file name)
    % removed handling of repeated entries for concatenated input data
       % for combining replicates, use 'combine_readTables.m'
    % auto-define xy limits on plots (script easily changed to set limits)
    % improved script annotation
    
%INPUTS
% <mark> is a vector
    % empty vector [] gives interactive Data Tips for coordinates on plot
    % vector of insertion sites can be used to mark scatterplot
%<sampleA>, <sampleB>, <filterF> = sampleName strings from fileName
    % Infiles are insertion site read count tables from INseq pipeline
    % filterF is a filtered set of insertion sites used to filter A and B coordinates

%OUTPUT
% Insertion sites in filterF absent from A and/or B, added with totRead=0
    % then +1 totRead added to ALL sites prior to plotting
% Scatterplot of insertion site read counts in A-vs-B (not auto saved)
    % user can define sites to be marked (red fill) in function call
    % OR explore plot using cursor and DataTips window
        % X and Y = read count for samples A and B at Site (see below)
        % i = index of Site (ignore)
        % Site = genome coordinate of Tn insertion site

%% Collect files and data
%input list of file names, held in structure
files = dir('INSEQ*processed*'); 
%reorder from UNIX to natural numeric sequence
files = natsortfiles({files.name});
%remove files that aren't read tables
idx = ~contains(files,'filter_cpm');
files = files(idx);

%% 
function [plotSampleName out_data] = import(sample,files)
%collect sample files (specific, so 'A1' doesn't retrieve 'A11', 'A12')
fileName = files(contains(files,strcat('.scarf_',sample,'.bowtiemap')));
fileName = cell2mat(fileName);
%Import data tables
tempName='tempFile.txt';
copyfile(fileName,tempName);
out_data = dlmread(tempName,'','B1..E(end)');
out_data = sortrows(out_data,1);
delete(tempName)
%Replace underscore in sampleName with dash for plotting
plotSampleName = sample;
und = strfind(sample,'_');
plotSampleName(und)='-';
end

%% 
[pSampleA, clrtA] = import(sampleA,files);
[pSampleB, clrtB] = import(sampleB,files);
[pSampleF, clrtF] = import(filterF,files);
filter = clrtF(:,1);

%Get Coordinates for A and B that are in F (filter)
%filter = clrtC(:,1);
tf = ismember(clrtA(:,1),filter,'rows');
ctA = clrtA(tf,1);
fclrtA = clrtA(tf,:);
tf = ismember(clrtB(:,1),filter,'rows');
ctB = clrtB(tf,1);
fclrtB = clrtB(tf,:);

%Get coordinates in F that are NOT in A and B
inFnotA = setdiff(filter,ctA);
inFnotB = setdiff(filter,ctB);

%Add back missing coordinates at totRead=0
fill = zeros(length(inFnotA),3);
all_A = [fclrtA; inFnotA fill];
all_A = sortrows(all_A,1);
fill = zeros(length(inFnotB),3);
all_B = [fclrtB; inFnotB fill];
all_B = sortrows(all_B,1);
%Add one totRead to all coordinates in both A and B
all_A(:,4) = all_A(:,4) +1;
all_B(:,4) = all_B(:,4) +1;

%get TotalReads for coordinates shared by A and B | in F
abc = intersect(all_A(:,1), all_B(:,1)); %coordinates in common
trA = all_A(:,4);
trB = all_B(:,4);
trAB = trA./trB;
out = [abc trA trB trAB];

%PLOTTING
f1=figure;
scatter(trA,trB, 72,'MarkerEdgeColor','blue','DisplayName','TNinsertions'); 
set(gca,'xscale','log','yscale','log');
%xlim([1 1e+07]);
%ylim([1 1e+07]);
set(gca, 'FontSize', 18);
name = strcat(pSampleA,'.vs.',pSampleB,'-CoordFrom:',pSampleF);
title(name,'FontSize', 16);
xlabel(pSampleA,'FontSize', 20);
ylabel(pSampleB,'FontSize', 20);

%Use 'mark' input to determine behavior:
if isempty(mark)==1
    %Show coordinates in Data Tips window
    coor = cellstr(num2str(abc));%must be cell array
    %interactive gene viewer
    dcm_obj = datacursormode(f1);
    set(dcm_obj,'UpdateFcn',{@showCoordinate,coor})
    datacursormode on
else
    %Add labels to specific coordinates defined by user input
    hold on
    [tf,mind] = ismember(mark,all_A);
    scatter(trA(mind),trB(mind),54,'filled','MarkerFaceColor','r',...
        'DisplayName','markedByUser');
    legend('Location','northwest');
    %position([0.177 0.749 0.163 0.129]);
    hold off
end

%figfile = strcat(sampleA,'_vs_',sampleB,'_V2filter_',filterF);
%savefig(f1,figfile);
end








