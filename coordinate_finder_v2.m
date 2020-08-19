function out=coordinate_finder_v2(LRout,totReadThresh,coordinates)
%last edit, July-04-2019
%REQUIRES natsortfiles.m

%OVERVIEW, version-2
    % Load read tables from current directory into a structure
        %can be used in Goodman 'results' directory, selects for readTables
    % Calculate scaling factor for proportional normalization to CPM
        % Omit coordinates w/ L/R > 1000, to remove phony read-suckers
    % Apply user defined threshold for total read count
    % Retrieve user defined coordinates and data
    % Apply scaling factor
%INPUTS
    % 'LRout' flag: 0==NO; 1==output L & R reads in data table (before Tot)
    % 'totReadThresh' is threshold for total reads, integer
    % 'coordinates' vector of TN insertion coordinates   
  % e.g. >> weirdo=coordinate_finder_v2(1,10,1459387)
%OUTPUT table (sorted by coordinate):
    % sample / coordinate /         / Tot-read / normalized Tot-read (CPM)
                   % / L-read / R-read /
                   
%% Consolidate sample names & read data into structure
%input list of file names, held in structure
files = dir('INSEQ*processed*'); 
%reorder from UNIX to natural numeric sequence
files = natsortfiles({files.name});
%remove files in Goodman 'results' directory that aren't read tables
idx=~contains(files,'filter_cpm');
files=files(idx);
%save sample name and data into structure
s = struct('sample',{},'data',[]);
tempName = 'tempFile.txt';
for i=1:length(files)
    infile = char(files(i));
    copyfile(infile, tempName);
    s(i).sampleID = infile(strfind(infile,'scarf_')+6:strfind(infile,'.bowtiemap')-1);
    s(i).data = dlmread(tempName,'',0,1);%column offset=1 to skip genomeName
end
delete(tempName)

%% Filter data, collect coordinates, build output
out={}; %create cell array for output
%Apply filterFile: keep coordinates of infile that are in filterFile
for i=1:size(s,2)
    clrt=s(i).data;
    
    %Filter & calculate scaling factor
    %remove rows with totRead < initial cutoff for total reads (=3)
    carCut = clrt(clrt(:,4) >= 3, :);
    %remove rows with L or R < cutoff (minLR)
    minLR=1; cutRat=1000;
    carCut_noz = carCut((carCut(:,2)>=minLR),:);
    carCut_noz = carCut_noz((carCut_noz(:,3)>=minLR),:);
    %calculate ratio of L-to-R, apply cutoff (cutRat)
    ratLR=carCut_noz(:,2)./carCut_noz(:,3);
    ccnzr = [carCut_noz ratLR];
    ccnzr = ccnzr((ccnzr(:,5)>=(1/cutRat)),:);
    ccnzr = ccnzr((ccnzr(:,5)<=(cutRat)),:);
    %calculate scaling factor
    tot=sum(ccnzr(:,4));
    norm2cpm=10^6/tot;
    
    %Find coordinates & generate output
    tf = ismember(clrt(:,1),coordinates);
    fclrt = clrt(tf,:);
    if LRout==1 %retrieve Land R read data
        hold=fclrt; totCol=4;
    else %just coordinate and totReads
        hold=fclrt(:,[1 4]); totCol=2;
    end
    %apply read threshold
    hold=hold(hold(:,totCol)>totReadThresh,:);
    %apply scaling factor, add CPM data
    cpm=hold(:,totCol) .* norm2cpm;
    hold=[hold cpm];
    %convert to cell array
    chold=cell(size(hold,1),1);
    chold(:,1)={s(i).sampleID};
    chold=[chold num2cell(hold)];
    out = [out;chold];
end
%sort output table by genome coordinate
out=sortrows(out,2);
end


