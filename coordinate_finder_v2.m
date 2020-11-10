function out=coordinate_finder_v2(LRout,totReadThresh,coordinates)
%last edit, Nov-11-2020 (annotation update)
    % requires "natsortfiles.m"
%e.g. >> myFavoriteSites = coordinate_finder_v2(1,10,[89587, 1459383])

%INPUTS
    % 'LRout' flag: 0==NO; 1==output L & R reads in data table (before Tot)
    % 'totReadThresh' is threshold for total reads, integer
    % 'coordinates' vector of TN insertion coordinates (# in sq.brackets)
%OUTPUT table (sorted by coordinate):
    % sample / coordinate / [L-read / R-read] / Tot-read / CPM

% version-2, additional from v1
    % Load read tables from current directory into a structure
        %can be used in Goodman 'results' directory, selects for readTables
    % Calculate scaling factor for proportional normalization to CPM
        % For CPM, apply L/R filters (minLR=1, ratioLR > 1000)
    % Retrieve user defined coordinates and data
        % Apply user defined threshold for total read count
        % Apply scaling factor
                   
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
for i=1:size(s,2)
    clrt=s(i).data;
    
    %Calculate scaling factor (after noise filter)
    %remove rows with totRead < initial cutoff for total reads (=3)
    carCut = clrt(clrt(:,4) >= 3, :);
    %remove rows with L or R < cutoff (minLR)
    minLR=1; cutRat=1000; %junk filter
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
    
    %Find coordinates & generate output (unfiltered)
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


