% July-2020
% For saving filtering list in Goodman format

% DEFINE workspace variable 'take_data'
    % nx4 array of TAsite/L/R/total counts
take_data = no8_newTAdata_omit8;
% CHOOSE passfile name
passfile='INSEQ_experiment.scarf_outside8-newTA-wNo8.bowtiemap_processed.txt_BEI_concat';
% CHOOSE chrID name    
chrID={'>BEI_concat'};

chrID(1:length(take_data))=chrID;
chrID=chrID';
passcell = [chrID num2cell(take_data)];
% write cell array to tab-delimted text
fileID = fopen(passfile,'w');
formatSpec = '%s\t %d\t %d\t %d\t %d\n';
[nrows,ncols] = size(passcell);
for row = 1:nrows
    fprintf(fileID,formatSpec,passcell{row,:});
end
fclose(fileID);