
%INPUTs: Excel spreadsheet (see content & formatting detail below)
    %create cell array  >>input={};
    %paste in Excel data (4 columns: site, plate, row, column)
%OUTPUT: 
    %histogram of TN#
    %cell arrays with information
        %'oneTN_out' & 'multiTN_out'
        %coordinate; #TN per well; other coordinates in well; well-ID
    %sort by coordinate# and save to tab-delimited Table?

%% Parse input, build binary for coordinate-well map

%Excel Sudoku tables have mix of numbers and letters, separated by commas
    %Matlab doesn't like reading these; numbers get truncated
%in Excel, replace commas with 'x', then import as cell array, 'input'
allSites=cell2mat(input(:,1));
sudin=input(:,2:4);
%Convert to cell array of strings, each element separated
sudout=cell(size(sudin));%for generating logicals
sudout_form=cell(size(sudin));%formatted for output to Excel
%sudlog=cell(size(sudin));
for i=1:size(sudin,2)
    for j=1:size(sudin,1)
        if isnumeric(sudin{j,i})
            sudout(j,i) = {num2str(sudin{j,i})};
            sudout_form(j,i) = sudout(j,i);
        else
            sudout(j,i) = {strsplit(sudin{j,i},'x')};
            sudout_form(j,i) = {strjoin(sudout{j,i},',')};
        end
    end
end

%Reference lists of elements
plates={'23','H','J','17','5','L','8','G','50','7','24','13','E',...
    '21','9','26','4','K','51'};
rows={'A','B','C','D','E','F','G','H'};
cols={'1','2','3','4','5','6','7','8','9','10','11','12'};

%Generate logicals for each element
sudlog=cell(size(sudin));
for i=1:size(sudin,2)
    if i==1, refer = plates;
        elseif i==2, refer = rows;
        elseif i==3, refer = cols;
    end
    for j=1:size(sudin,1)
        sudlog{j,i}=ismember(refer,sudout{j,i});
    end
end

%Using logicals to extract Tn# information
%combine into single logical array of elements: 19+8+12 = 39 columns
all_sudlog=cell2mat(sudlog);
%
[C,ia,ic]=unique(all_sudlog,'rows','stable'); %all_sudlog=C(ic,:)
%histogram of ic values (indices of unique and repeating elements)
[A,B]=hist(ic,unique(ic));
iM = ~ismember(ic,B(A>1));%logical: one-TN=1, multi-TN=0

%Summary stats
  %distribution of #Tn/well-ID for input
    [cnt,nTN]=hist(A,unique(A));
 
%% Collect information for one-TN entries
%Use logical to select one-TN data
oneTN_sites=num2cell(allSites(iM));%cell
fill1=num2cell(ones(length(oneTN_sites),1));
fill2=cell(length(oneTN_sites),1);
oneTN_well=sudout_form(iM,:);%cell
%combine for output
oneTN_out=[oneTN_sites,fill1,fill2,oneTN_well];

%% Collect information for multi-TN entries
%'ic' contains index for first instance of coordinate
  %retrieve all instances and extract:
    %coordinate
    %number of coordinates that share well-ID (#Tn)
    %vector of shared coordinates
    %shared well-ID

%first index of replicate entries
fim=find(hist(ic,unique(ic))>1);

mTN_sites={};
mTN_sud={};
for i=1:length(fim) %collect replicates, sequentially
    aim = find(ic==fim(i)); %all indices for multi-TN (i)
    hold_sites=cell(length(aim),3);
    hold_sud=cell(length(aim),3);
    for j=1:length(aim)%
        hold_sites{j,1}=allSites(aim(j));%coordinate for multi-TN
        hold_sites{j,2}=length(aim);%number that map to well
        other= aim(aim(j)~=aim);
        temp=cell(length(other),1);
        for k= 1:length(other)
            %other coordinates that map to well
            temp(k) = {num2str(allSites(other(k)))};
        end
        if length(temp)>1
            hold_sites{j,3}={strjoin(temp(:),',')};
        else
            hold_sites{j,3}=temp;
        end
        hold_sud(j,:)=sudout_form(aim(1),:);%well-ID for multi-TN
    end
    mTN_sites=[mTN_sites;hold_sites];
    mTN_sud=[mTN_sud;hold_sud];
end
mTN_out=[mTN_sites,mTN_sud];
    