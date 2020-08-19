function filter_readTable_byTAsites(sampleName,TAsites)
%last edit: July-19-2020

% Filters each using user-defined coordinates for TAsites and saves
   % filtered read tables to new directory "TA-filtered_readTables"
    
%INPUT   
% User-defines action on single file or directory
   %for single file, <sampleName> = string of sample name in file name
      %filter_readTable_byTAsites('Input1','Akk_TAsites.txt')
   %for directory, <sampleName> = any number
      %filter_readTable_byTAsites(1,'Akk_TAsites.txt')
% Goodman-format read tables
% <TAsites> is string for .txt list of genome coordinates (Artist format)

%% Directory and file handling

%Import all TA sites from TAsites.txt file
TAcoor=dlmread(TAsites,'',0,1); %1 skips first column

%Process input files
files = dir('INSEQ*processed*'); %input list of file names in structure
%Reorder from UNIX numeric sequence to standard
files = natsortfiles({files.name});
%remove files in Goodman 'results' directory that aren't read tables
idx = ~contains(files,'filter_cpm');
files = files(idx);
%restrict to single file based on user input
if isa(sampleName,'char')==1
    idx=contains(files,strcat('scarf_',sampleName,'.bowtiemap'));
    files=files(idx);
end

% Write to new directory, check for overwrite
TAdir='TA-filtered_readTables';
if exist(TAdir,'dir')==0
    mkdir(TAdir);
else
    disp('"TA-filtered_readTables" directory already exists');
    m=input('Do you wish to continue, y/n:','s');
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
    %Remove non-TA sites from data file
    indata=indata(ismember(indata(:,1),TAcoor),:);
    
    %Save output
    %read character column from infile and store as cell array for output
    T = readtable(tempName, 'ReadVariableNames',false, 'Delimiter', 'tab');
    chrID = table2array(T(:,1)); %resize before output
    chrID = chrID(1:size(indata,1));
    passcell = [chrID num2cell(indata)];
    passfile = strcat(TAdir,filesep,infile);
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

