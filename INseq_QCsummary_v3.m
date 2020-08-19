function INseq_QCsummary_v3(plot,filterName)
%last edit, July-6-2019

% e.g. >> INseq_QCsummary_v3(1,'ArrayPoolSelect')
    
%Acquires metrics from read tables to evaluate quality/content
%changes for v3: build QC table for all sample files in current directory 
    % curDir contents = multiple sample files and single filter file
    % filter file <SAMPLE> name must be unique (=filterName input string)
    % removed metrics output in v2:
        % L/R min & max
        % <outliers> for 'moreThan3' filtered data
%from v2 description:
    % only plots filter C data (filterFile)
    % L/R metrics come from log space
    % top15 L/R outliers extracted with metric converted to >1
        %coordinate and L/R or R/L value (whichever is >=1)

%Overview, data handling
    %Filter-A: none
    %Filter-B: select coordinates with totReads >3 (Goodman default)
    %Filter-C: select insertion coordinates of "filterFile"
    %Filter-D: Filter-C data, then apply L/R threshold
    %note: for L/R, add pseudocount to L or R=0 to allow L/R calc. for all
        %apply after filtering & extracting other read metrics
        %done by replacing L,R=0 value with 1  
%INPUT
    % <plot> is flag: 0, no; 1, yes
    % <filterName> is a string unique to the file name used for filtering
        % e.g. 'BigArrayPool' (not the whole file name)
    % filtering file should be pre-filtered by user (e.g. reads, knowns)
    % data table in current directory should have standard naming:
       %"INSEQ_experiment.scarf_<SAMPLE>.bowtiemap_processed.txt_<GENOME>"
%OUTPUT is table of values for raw and filtered read tables
    % filter name
    % total #reads
    % fraction total raw reads (raw = unfiltered total reads)
    % #TN insertions (coordinates)
    % mean #reads per coordinate
    % median #reads per coordinate
    % L/R stats (add pseudocount to capture all coord., mark red on plot):
        %mean, median, variance 
        %outliers: fiftteen most extreme coordinate, then L/R,R/L value 
%PLOTS (optional, user inout; Filter-C only)
    % scatterplot of #totReads -vs- L/R for all coordinates (v3plot)
    % pseudocount=1 added for L or R=0,to capture L/R for these coor.
        % coordinates with +pseudo marked on plot in red
    % figure saved by default
    
%% Consolidate sample data in a structure
files = dir('INSEQ*processed*'); %input list of file names in structure
%Reorder from UNIX numeric sequence to standard
files = natsortfiles({files.name});
%remove files in Goodman 'results' directory that aren't read tables
idx=~contains(files,'filter_cpm');
files=files(idx);
%Keep coordinates from filtering file in a vector outside stucture
idx=contains(files,filterName);
filterFile=char(files(idx));
files=files(~idx);
tempName = 'tempFile.txt';
copyfile(filterFile, tempName);
clrtF = dlmread(tempName,'',0,1);
%clrtF = sortrows(clrtF,1);
filter = sort(clrtF(:,1));
%Load structure
s = struct('sample',{},'data',[]);
for i=1:length(files)
    infile = char(files(i));
    copyfile(infile, tempName);
    sampleName = infile(strfind(infile,'scarf_')+6:strfind(infile,'.bowtiemap')-1);
    s(i).sample=cellstr(sampleName);
    s(i).data = dlmread(tempName,'',0,1);%skip col.1 genomeName  
end
delete(tempName) 

%% Retrieve data, apply filters, extract metrics, hold in structure
for i=1:size(s,2) 
    clrt=s(i).data;
    sampleName=s(i).sample;
    %FILTER-B: Goodman-style, totReads > 3
    totReads = 3;% Note:Goodman script uses >3 total reads filter (not >=3)
    clrtB = clrt((clrt(:,4) > totReads), :);
    %FILTER-C: keep coordinates of infile that are in the filtering file
    tf = ismember(clrt(:,1),filter);
    clrtC = clrt(tf,:);
    %FILTER-D: Filter-C plus cutoff for L/R bias
    clrtD=clrtC; %start with coordinate filtered data
    %Add temporary PSEUDO to LorR=0 to calculate L/R for filtering 
    Li=find(~clrtD(:,2)); %indices with L=0;
    Ri=find(~clrtD(:,3)); %indices with R=0;
    clrtD(Li,2)=1; clrtD(Ri,3)=1; %replace 0 with 1
    LRrat = clrtD(:,2) ./ clrtD(:,3);
    clrtD=clrtC(LRrat>0.01 & LRrat<100,:);%final data lacks pseudocounts
    %Hold filtered data in structure
    s(i).filtTotRead=clrtB;
    s(i).filtFile=clrtC;
    s(i).filtFileLRcut=clrtD;
    
    %Collect metrics for different filters, consolidate output
    total=sum(clrt(:,4));%all reads
    r1 = getNum(clrt,total); %INTERNAL FUNCTION (see below)
    r2 = getNum(clrtB,total);
    r3 = getNum(clrtC,total);
    r4 = getNum(clrtD,total);
    %build OUT array: combine output rows; add filter names
    nout=[r1;r2;r3;r4];
    nout = num2cell(nout);
    fout = cell(4,1);
    fout{1}=sampleName; fout{2}='moreThan3';fout{3}=filterName;
    fout{4}=strcat(filterName,'-LRcut');
    fout(:,2:6)=nout;   

    %Collect L/R bias information; generate scatter plot using file filter
    [lovrB,~] = getLR(0,clrtB,sampleName,'moreThan3');
    [lovrC,lierC] = getLR(plot,clrtC,sampleName,filterName);
    Dname=strcat(filterName,'-cutLR100');
    [lovrD,~] = getLR(0,clrtD,sampleName,Dname);
    %output cell array
    fill=cell(1,length(lovrB));
    LRcell=[fill;lovrB;lovrC;lovrD];
    fout = [fout LRcell];
    fout=fout';
    outliers=lierC;
    
    %reshape for easier output
    rout=[fout(1:6,1);fout(1:9,2);fout(1:9,3);fout(1:9,4);];
    
    %hold in structure
    s(i).metrics=rout; %cell array
    s(i).outliers=num2cell(outliers); %cell array
    
end

%% Consolidate output metrics & generate outfile

outTable=[s.metrics];
label(1,1:size(outTable,2))={'outliers, insertion site'};
outTable=[outTable;label];
outies=[s.outliers];
label(1,1:size(outTable,2))={'outliers,  L/R or R/L ratio'};
outTable=[outTable; outies(1:15,:); label; outies(16:30,:)];

%collect sample names from cell array for column headers
%variableNames = outTable(1,:); NEED to reformat for 'writetable'
%outTable=outTable(2:end,:);
outfile = 'QC_v3_table';
Tout = cell2table(outTable);
writetable(Tout,outfile,'FileType','text','Delimiter','tab');


%FUNCTIONS:
%% Extract metrics from all filtrates
function out=getNum(inData,total)
%total=sum(clrt(:,4));%all reads
totX=sum(inData(:,4));%all reads after filterX
fractTot=totX/total;%percent of unfiltered reads in X
tns=size(inData,1);%number of unique insertion sites
meanRX=mean(inData(:,4));
medRX=median(inData(:,4));
%generate output vector
out=[totX fractTot tns meanRX medRX];
end

%% Extract L/R bias information
function [out,outliers]=getLR(plot,inData,sampleName,filterName)
%<plot> input varaible is flag for generating plot: 0=no, 1=yes
%ADD PSEUDO to capture completeley lopsided data (e.g. L=0, R=100)
  %hold on to +PSEUDO coordinates, for marking on plots, etc...
zLi=find(~inData(:,2)); %indices with L=0;
zRi=find(~inData(:,3)); %indices with R=0;
inData(zLi,2)=1; inData(zRi,3)=1; %replace 0 with 1
%calculate ratio of L-to-R (no cutoff)
LovR = inData(:,2) ./ inData(:,3);

%Convert to all >=1 (if L/R <1, take R/L):
LRdata= [inData(:,1), LovR];%coordinate & L/R
LRdata=sortrows(LRdata,2);
flip=find(LRdata(:,2)<1);
LRg1=LRdata;
LRg1(flip,2)=1./(LRg1(flip,2));%converted to >=1
%Retrieve coordinate and >1 metric (L/R or R/L) for top10 L/R bias
[LRg1,ind]=sortrows(LRg1,2,'descend');
top15=LRg1(1:15,:);
%top10=LRdata(ind(1:10),:);%retrieve L/R
outliers=vertcat(top15(:,1),top15(:,2));
%maybe: change to cell array, add filterName header

%Get all statistics from log10 data:
logLR=log10(LovR);
meanLovR=mean(logLR);
medLovR=median(logLR);
varLovR=var(logLR);
out=[meanLovR medLovR varLovR];
out=num2cell(out);
%[l,h]=bounds(logLR);           SUPRESS min & max output
%bndLovR=num2cell([l h]);
%out = [out bndLovR];

%Plot total reads vs. ratio of LtoR 
if plot==1
    figure;
    scatter(inData(:,4),LovR);
    set(gca,'Xscale','log','Yscale','log');
    set(gca, 'FontSize', 12)
    name=strcat(sampleName,'--',filterName);
    title(name);
    xlabel('total mapped reads (per coordinate)');
    ylabel('ratio of L-to-R mapped reads (per coordinate)');
    %mark coordinates that had L or R=0, added pseudocount=1
    mark=[zLi;zRi];
    hold on
    scatter(inData(mark,4),LovR(mark),28,'filled','MarkerFaceColor','r',...
            'DisplayName','pseudo_added');
    hold off
    %save file?
    figName=cell2mat(strcat(sampleName,'_filt-',filterName,'_QCv2.fig'));
    savefig(figName);
    close
    
end
end

end


