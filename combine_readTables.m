function combine_readTables(mode,samples)
%last edit: July-17-2020

%Run from within current directory
%Consolidates data from multiple read table files defined by user
	%Improved how duplicates are handled for efficiency - MUCH faster!

%INPUT
    %<mode>, string 'sum' or 'mean'
    %<sampleIDs>, cell array of sample names for read tables (Goodman output)
        %e.g. {'sample1','sample2'}
        %files contents = Genome//Coordinate//#Reads(L/R/Tot)
%OUTPUT
  %Consolidated data table, same format as inputs
  %New name is concatenation of <mode> and sample names

%% Collect files and data
%input list of file names, held in structure
files = dir('INSEQ*processed*'); 
%reorder from UNIX to natural numeric sequence
files = natsortfiles({files.name});
%collect sample files
 %more specific, so 'A1' doesn't retrieve 'A11', 'A12' (July-17-2020)
idx=contains(files,strcat('.scarf_',samples,'.bowtiemap')); 
files=files(idx);
%remove matching files that aren't read tables
idx=~contains(files,'filter_cpm');
files=files(idx);

%Check before combining
if length(samples)~=length(files)
    disp('SCRIPT TERMINATED!');
    disp('Number of samples defined by user does not match #of files found');
    return
end

%loop through all sample files
tempName='tempFile.txt';
s = struct('data',[]);
for i=1:length(files)
    infile = char(files(i));
    copyfile(infile, tempName);
    %read numeric columns from infile, sort by insertion coordinate
    clrt = dlmread(tempName,'','B1..E(end)');
    clrt = sortrows(clrt,1);
    s(i).data = clrt;
    %
end
%read character column from infile and store as cell array for output
T = readtable(tempName, 'ReadVariableNames',false, 'Delimiter', 'tab');
chrID = table2array(T(1,1)); %resize before output
delete(tempName);

%% Aggregate data, find unique sites, apply <mode> function
all_data = cat(1,s.data);
%Handle multiple entries for chr.coordinate
  %sum or mean of multiples based on input flag <mode>
[xnew,~,idx] = unique(all_data(:,1));
if strcmp(mode,'mean')==1
    Lnew = accumarray(idx(:),all_data(:,2),[],@mean);
    Rnew = accumarray(idx(:),all_data(:,3),[],@mean);
    Tnew = accumarray(idx(:),all_data(:,4),[],@mean);
elseif strcmp(mode,'sum')==1
    Lnew = accumarray(idx(:),all_data(:,2));
    Rnew = accumarray(idx(:),all_data(:,3));
    Tnew = accumarray(idx(:),all_data(:,4));
else
    disp('define MODE for combining reads from multiple files');
    return
end
uData=[xnew Lnew Rnew Tnew];
uclrt = round(uData);

%OUTPUT consolidated data table in Goodman format
% Save output read table
newName = strcat(mode,'-',strjoin(samples,'-'));
  %without separators: newName = strcat(mode,cell2mat(samples));
chrName=cell2mat(chrID);
chrName=chrName(2:end);
passfile = strcat('INSEQ_experiment.scarf_',newName,'.bowtiemap_processed.txt_',chrName);
%resize cell array, incorporate filtered data
chrID(1:length(uclrt)) = chrID;
chrID = chrID';
passcell = [chrID num2cell(uclrt(:,1:4))];
%write cell array to tab-delimted text
fileID = fopen(passfile,'w');
formatSpec = '%s\t %d\t %d\t %d\t %d\n';
[nrows,ncols] = size(passcell);
for row = 1:nrows
    fprintf(fileID,formatSpec,passcell{row,:});
end
fclose(fileID);
end

