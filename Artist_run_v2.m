function Artist_run_v2(genome,inputName,sampleName,filterName)
%last edit: July-03-2020     %works

% v2 changes:
    % 'my' TAsites-filtered output only
    % don't round mean of controlsims (myNormInput)
    % [re]scale myNormInput to mySample (depth variation after simulation)
    
% Consolidates ARTIST workflow
% Generate normalized Input (mean of ControlSims)

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
% OUTPUT
    % .mat file named by inputName_sampleName_filteName
    % contains workspace variables for analysis and MWU stats
        % also TAhits: Tn insertions per locus
        % and NormInput: mean of simulated/normalized counts for each site


% Call set-up script
Artist_setup(genome,inputName,sampleName,filterName);
saveName = strcat(inputName,'_',sampleName,'_',filterName,'.mat');
load(saveName);

% Generate uniquenames and indices
[myUniqueNames,myUniqueIndices]=getuniquetanames(myTAid);
% Run simulations for normalization
[myControlSims]=simulateEqualSaturation_v3(myInput,mySample,100);
% Run Mann-Whitney test on annotated loci for each simulation
[myMWUboots] = runmwuallboots(myControlSims,mySample,myUniqueNames,myUniqueIndices);
% Generate stats
[myMWUstats]=MWUsummarystats(myMWUboots,0.05,0.9,myUniqueIndices,myControlSims,mySample);
% Generate <TAhits>
mui=myUniqueIndices;
mus=mui(:,2)-mui(:,1);
myTAhits=mus+1;
% Generate <NormInput> varirable
myNormInput=mean(myControlSims,2);
% [re]match depth to sample (stochastic changes from simulation)
myNormInput = myNormInput .* sum(mySample)/sum(myNormInput);

% Save all workspace variables
saveName = strcat('MWUdone-',saveName);
save(saveName,'TA*','Input','Sample','my*','genome',...
    'inputName','sampleName','filterName');

% Save data tables
passfile=strcat('MWUstats-',inputName,'_',sampleName,'_',filterName,'.txt');
%resize cell array, incorporate filtered data
out = [myTAhits myMWUstats];
passcell = [myUniqueNames num2cell(out)];
%write cell array to tab-delimted text
fileID = fopen(passfile,'w');
formatSpec = '%s\t %d\t %d\t %d\t %d\t %d\t %d\t %d\n';
[nrows,ncols] = size(passcell);
for row = 1:nrows
    fprintf(fileID,formatSpec,passcell{row,:});
end
fclose(fileID);

