function filter_readTable_byFile(sampleName,filterName)
%last edit: August-30-2020
%>> filter_readTable_byFile('GF7d-cecum1','Inputs-tot20')
%>> filter_readTable_byFile(1,'Inputs-tot20')

% Filters read table(s) with user-defined coordinates from file and
  % saves filtered read table(s) to new sub-directory
    
%INPUT   
  % User-defines action on single file or directory
    %for single file, <sampleName> = string of sample name in file name
    %for directory, <sampleName> = any number
  % Goodman-formatted read tables
%OUTPUT
  % sub-directory "<filterName>-filtered_readTables"
    % containing filtered read tables for all data in current directory 

%% Directory and file handling
%Process input files
files = dir('INSEQ*processed*'); %input list of file names into structure
%Reorder from UNIX numeric sequence to standard
files = natsortfiles({files.name});
%remove files in Goodman 'results' directory that aren't read count tables
idx = ~contains(files,'filter_cpm');
files = files(idx);
%retrieve filterName file and Tn sites (filter)
idx = contains(files,filterName);
filterFile = char(files(idx));
tempName = 'tempFile.txt';
copyfile(filterFile,tempName);
%remove filterFile from file list
files = files(~idx);
%import filterSites list
clrtF = dlmread(tempName,'',0,1);
filter = sort(clrtF(:,1));

%restrict action to single file based on user input
if isa(sampleName,'char') == 1
    idx = contains(files,strcat('scarf_',sampleName,'.bowtiemap'));
    files = files(idx);
end

% Write to new directory, check for overwrite
filterDir = strcat(filterName,'_filtered_readTables');
if exist(filterDir,'dir') ==0
    mkdir(filterDir);
else
    msg = strcat('"',filterDir,'" directory already exists');
    disp(msg);
    m = input('Do you wish to continue and overwrite, y/n:','s');
    if m=='n'
        return
    end
end

%% Filtering and output  
tempName = 'tempFile.txt';
for i=1:length(files)
    %Extract data tables
    infile = char(files(i));
    copyfile(infile,tempName);
    indata = dlmread(tempName,'',0,1);%column offset=1 to skip genomeName
    indata = sortrows(indata,1);
    %Remove sites (rows) from data file that are not in filterFile
    outdata=indata(ismember(indata(:,1),filter),:);
    
    %Save output
    %read character column from infile and store as cell array for output
    T = readtable(tempName, 'ReadVariableNames',false, 'Delimiter', 'tab');
    chrID = table2array(T(:,1)); %resize before output
    chrID = chrID(1:size(outdata,1));
    passcell = [chrID num2cell(outdata)];
    passfile = strcat(filterDir,filesep,infile);
    %write cell array to tab-delimted text
    fileID = fopen(passfile,'w');
    formatSpec = '%s\t %d\t %d\t %d\t %d\n';
    [nrows,ncols] = size(passcell);
    for row = 1:nrows
        fprintf(fileID,formatSpec,passcell{row,:});
    end
    fclose(fileID);
end
delete(tempName);
end
