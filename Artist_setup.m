function Artist_setup(genome,inputName,sampleName,filterName)
%last edit: April-22-2020
%works

% Parse all inputs and generate .mat file containing all workspace
    % variables needed for conditional gene analysis with ARTIST 

% INPUTs
    % 'Genome', prefix for Genome_TAsites and Genome_TAsitesID .txt files
    % 'inputName', filname: sample name for read table from Goodman pipeline
        % coordinates in column-2, total reads in column-5
    % 'sampleName', filname: sample name for read table from Goodman pipeline
        % coordinates in column-2, total reads in column-5
    % OPTIONAL: 'filterName', sample name for read table from Goodman pipeline
        % used to filter sample data prior to filling missing TA sites
        % leave out of command to use input data table as is
        
% >> Artist_setup('Akk','GF7dinoc','CecumT1F1','BAPfeb19')

% OUTPUT
    % .mat file named by inputName_sampleName_filteName
    % contains workspace variables requirerd for ARTIST
        % <TAsites>
        % <TAid>
        % <Input>
        % <Sample>
        % AND additional versions of the same with the prefix 'my' that are
          % filtered for coordinates found in the filterName file, if given
        
%Import all TA sites from Genome_TAsites.txt file
TAsites_file=strcat(genome,'_TASites.txt');
T=readtable(TAsites_file);
TAsites=table2array(T(:,2));

%Import all TAid from Genome_TAsitesID.txt file
TAid_file=strcat(genome,'_TAsitesID.txt');
T=readtable(TAid_file);
TAid=table2cell(T(:,2));

%Access files
files = dir('INSEQ*processed*'); %input list of file names in structure
%Reorder from UNIX numeric sequence to standard
files = natsortfiles({files.name});
%remove files in Goodman 'results' directory that aren't read tables
idx=~contains(files,'filter_cpm');
files=files(idx);
    
[Input, myInput] = process_reads(inputName);
[Sample, mySample] = process_reads(sampleName);

%% extract, filter and fill data
function [out, myOut] = process_reads(name)
%Load data tables from input file
infile=char(files(contains(files,name)));
%Extract data tables
tempName = 'tempFile.txt';
copyfile(infile, tempName);
data = dlmread(tempName,'',0,1);%column offset=1 to skip genomeName
data = sortrows(data,1);
data = data(:,[1 4]);%keep only coordinate and total reads

%Remove non-TA sites from data file
data=data(ismember(data(:,1),TAsites),:);

%If absent from TAsites, add coordinate to data with read count =0
    %requires 'TAsite' and 'data' (coordinates & counts)
inTAnotData = setdiff(TAsites,data(:,1));
fill = zeros(length(inTAnotData),1);
data = [data; inTAnotData fill];
out = sortrows(data,1);
out = out(:,2);%reads for all TA sites

%Use coordinates in filterName file to limit TA sites
if exist('filterName','var')
    %Keep coordinates from filtering file in a vector outside stucture
    idx=contains(files,filterName);
    filterFile=char(files(idx));
    copyfile(filterFile, tempName);
    clrtF = dlmread(tempName,'',0,1);
    filter = sort(clrtF(:,1));
    %Filter data based on filterFile
    myOut = data(ismember(data(:,1),filter),:);
    myOut = sortrows(myOut,1);
    myOut = myOut(:,2); %reads for selected TA sites
else myOut=[];
end
delete(tempName);
end

%% Apply filters to TA variables, save workspace variables for ARTIST
%Use coordinates in filterName file to limit TA sites
if exist('filterName','var')
    %Keep coordinates from filtering file in a vector outside stucture
    idx=contains(files,filterName);
    filterFile=char(files(idx));
    holdName = 'holdFile.txt';
    copyfile(filterFile, holdName);
    clrtF = dlmread(holdName,'',0,1);
    delete(holdName);
    filter = sort(clrtF(:,1));
    %Filter based on filterFile
    mySites = ismember(TAsites,filter);
    myTAsites = TAsites(mySites);
    myTAid = TAid(mySites);
else filterName='noFilter';
end

saveName = strcat(inputName,'_',sampleName,'_',filterName,'.mat');
save(saveName,'TA*','Input','Sample','my*','genome',...
    'inputName','sampleName','filterName');

end

