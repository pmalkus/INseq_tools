function out = sudoku_plate_compare_v1(meanReplicatesLimit,meanReplicatesHist)
%Last edit, Feb-26-2019 (annotations updated Dec-09-2020)
%e.g. >> out = sudoku_plate_compare_v1(3,2)

%Run from within current directory
%Depends on proper naming convention of INPUT files!!! e.g.:
%'INSEQ_experiment.scarf_Plate9.bowtiemap_processed.txt_Akkermansia_BAA835_rv'

%INPUT:
  % User-filtered read tables of Tn-insertion sites for each 96-well Plate
     % Important to establish quality criteria for 'real' Tn-insertions
  % <meanReplicatesLimit>, scalar
     % Threshold to stop running optimization routine
  % <meanReplicatesHist>, scalar
     % Value of <meanReplicates> (see below) used to indentify a specific
     % plate set for creating a histogram of <replicates>
     
%OUTPUT:
  % Sequence of Plates optimized for low clonal redundancy
  % <out> is a structure, each element holding information on a plate set
    % e.g. out(4) is the set of 4 plates w/ lowest # of replicates
    % out.plates holds Plate IDs from input file names
    % out.data holds Tn-insertion sites and #replicates (intra & inter -plate)
  % Fig-1: line plots showing trade-off between coverage & redundancy
    % number of unique Tn-site VS number of plates in set
    % number of replicates (mean) VS number of plates in set
  % Fig-2: histogram of <replicates> for plate set = <meanReplicatesHist>
  % Fig-3: count of <replicates> = 1, 2, 3 and 4 VS plate set size
  
%PROCESS:
%Consolidate data for each Plate input table into structure
%Generate 'replicates' varible for each Tn-insertion coordinate
    %first, WITHIN plates:
        %find # unique Tn-insertion per plate (X)
        %randomly assign remainder (96-X) with replacement
    %then, BETWEEN plates, add replcate number based on shared coordinate

%Optimization routine:
%  1. Find best plate based on lowest mean(replicates) == (Plate-S1)
%  2. Sample all remaining plates to find best plate to combine with S1, 
    % to generate a plate set that minimizes mean(replicates)
%  3. Repeat until mean(replicates) for plate set = <meanReplicatesLimit> 
%Information from each round preserved in output structure

%% Consolidate plate data in a structure, generate 'replicates' variable
files = dir('INSEQ*'); %input list of file names in structure
%Reorder from UNIX numeric sequence to standard
files = natsortfiles({files.name}); %changed 190628
s = struct('plate',{},'data',[]);
for i=1:length(files)
    infile = char(files(i));        %changed 190628
    s(i).plate = infile(strfind(infile,'Plate'):strfind(infile,'.bowtie')-1);
    %Copy input file with a new name understood by 'readtable'
    tempName = strcat(infile(1:5),'.txt');
    copyfile(infile, tempName);
    out = dlmread(tempName,'','B1..B(end)');%limit to coordinates
    delete(tempName)
    %Generate 'replicates' variable to assign for each coordinate
    replicates = ones(length(out),1);
    %Sample from list of coordinates, w/replacement, to assign redundancy
    add = randsample(length(out),96-length(out),true);
    for j=1:length(add)
        replicates(add(j)) = replicates(add(j)) +1;
    end
    out = [out replicates];
    %add to structure
    s(i).data = out; %two column vectors: coordinate & replicates
end

%% Find best combination of plates to minimize value of mean(replicates)
%Set-up
choices = 1:size(s,2);%indices of s
itaken =[];%indices of s already selected for plate aggregate
%Save in structure; each index of 'suko' will have an additional plate
suko = struct('plates',{},'data',[]);
j=1; 
meanRep=0;
%limit = 2; %set limit for collecting plates in output structure
while meanRep < meanReplicatesLimit && length(itaken) < length(choices)
    iremain = setdiff(choices,itaken);%indices of s remaining to sample
    if j==1, base =[]; else base = suko(j-1).data; end
    %Find next best plate
    score=[];%place to hold mean(replicates) metric for each test case
    for k=1:length(iremain)
        itest = iremain(k);
        test = s(itest).data;%retrieve data for plate to test
        btest = [base; test];%combine with plates already picked
        [uC,~,idC] = unique(btest(:,1));%find shared coordinates
        uR = accumarray(idC(:),btest(:,2));%for shared, sum 'replicates' value
        score(k) = mean(uR);
    end
    %hold index of iremain used for aggreagate with best score (val)
    [meanRep,idx] = min(score);
    itaken = [itaken iremain(idx)];
    %add chosen plate to suko structure
    take = cat(1,s(itaken).data);
    [uC,~,idC] = unique(take(:,1));%find shared coordinates
    uR = accumarray(idC(:),take(:,2));
    suko(j).data = [uC uR];
    suko(j).plates = cat(1,{s(itaken).plate});
    %if j==1, suko(j).plates = s(iremain(idx)).plate;
    %else suko(j).plates = [suko(j-1).plates; {s(iremain(idx)).plate}];%cell!!!!
    %end
    j=j+1;
end
out = suko;
%% Plotting, etc...
%Plot mean(replicates) -vs- platesTaken
%   co-plot unique insertion sites -vs- platesTaken
%Build data vectors
meanRep=[];nSites=[]; 
isuk = length(suko);
for i=1:isuk
    meanRep(i) = mean(suko(i).data(:,2));
    nSites(i) = length(suko(i).data(:,2));
end
f1 =figure;
hold on
yyaxis left
plot(1:isuk,meanRep);
xlabel('number of plates included');
ylabel('mean number of replicates');

yyaxis right
plot(1:isuk,nSites);
ylabel('number of unique insertion sites')
hold off

% histogram of 'replicates' at mean(replicates)=meanReplicatesHist (mrh)
f2 = figure;
%find the index in suko with mean(replicates) closest to mrh
mrh = meanReplicatesHist;
[~,idR] = min(abs(meanRep-mrh));%beware multiple at same value, rare
reps = suko(idR).data(:,2);
hist(reps,max(reps));
xlabel('number of replicates');
ylabel('count');
xticks(1:max(reps));

%
%plot (replicates=1_counts)/(replicates=2_counts)  -vs-  plates included
f3=figure;
wuns=zeros(1,isuk);twos=wuns;threes=wuns;fours=wuns;
for i=1:isuk
    reps2 = suko(i).data(:,2);
    repHist = hist(reps2,1:15);
    wuns(i) = repHist(1);
    twos(i) = repHist(2);
    threes(i) = repHist(3);
    fours(i) = repHist(4);
end 
hold on
plot(1:isuk,wuns,'DisplayName','Single');
plot(1:isuk,twos,'DisplayName','Double');
plot(1:isuk,threes,'DisplayName','Triple');
plot(1:isuk,fours,'DisplayName','Quadruple');
legend('Location','northwest');
hold off
%
end




