function out=Sudoku_colRow_v1
%last edit, June-24-2019

%Overview
    %Loads files from curDir and extract column/row identifiers from names
    %Filter data table using insertion coordinates of "filterFile"
    %Take coordinates with top read counts
        %for columns: 8 rows X 19 plates = 152 coordinates (estimate)
        %for columns: 8 rows X 19 plates = 228 coordinates (estimate)
    %Build cell array: 21 columns x N rows; N= coordinates in filterFile
        %col.1 = coordinates from filterFile (pre-screened Plate data)
        %col.2-13 = column data
        %col.14-21 = row data
%INPUTS
    %Infile must be named correctly
        %Following 'Sud': Column number; Row letter; "Filter" for filter



%% Consolidate sample data in a structure
files = dir('INSEQ*'); %input list of file names in structure
%Reorder from UNIX numeric sequence to standard
files = natsortfiles({files.name});
s = struct('sample',{},'data',[]);
tempName = 'tempFile.txt';
for i=1:length(files)
    infile = char(files(i));
    copyfile(infile, tempName);
    sampleName = infile(strfind(infile,'Sud')+3:strfind(infile,'.bowtiemap')-1);
    if contains(sampleName,'ilter')
        filter = dlmread(tempName,'','B1..B(end)');%just coordinates
        filter = sort(filter);
    else
        %structure grows with each iteration
        s(size(s,2)+1).sample=cellstr(sampleName);
        s(size(s,2)).data = dlmread(tempName,'',0,1);
    end
end
delete(tempName)
if size(s)~=20
    disp('something is wrong. must have 20 data files');
end 
%% Filter data, collect coordinates, build output

%create cell array for output
assign=cell(length(filter),20);
%Apply filterFile: keep coordinates of infile that are in filterFile
for i=1:20 %has to be 20!
    clrt=s(i).data;
    tf = ismember(clrt(:,1),filter);
    fclrt = clrt(tf,:);
    out = sortrows(fclrt,4,'descend');%sort filt. coor. by most total reads
    if i<13
        top = out(1:152,1);
    else
        top = out(1:228,1);
    end
    %assign sample name (col/row) to matching coordinates
    [~,match,~]=intersect(filter,top);
    assign(match,i)=s(i).sample;
end
%Reorder from UNIX numeric sequence to standard
%assign=permute(assign(

%Collapse cell array into column-assignment and row-assignment columns
colRow=cell(length(filter),2);
for i=1:length(filter)
    ci=assign(i,1:12);
    temp=ci(~cellfun('isempty',ci));
    colRow(i,1)=cellstr(strjoin(temp,','));
    %join row info
    ri=assign(i,13:20);
    temp=ri(~cellfun('isempty',ri));
    colRow(i,2)=cellstr(strjoin(temp,','));
end   
%Add filter coordinates to output cell array
out=num2cell(filter);
out=[out, colRow];

end


