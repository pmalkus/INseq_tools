function INseq_coverage_cluster

%For oBAA-835
%linear bar graph: coordinate vs. total reads, passing filters
figure;
edges = [1:10000:2700000];
hist(CTcut(:,1),edges,'FaceColor',[0.8 0.2 0,2]); 
xlim([-40000 2700000]);
%xticks([0:3e+05:27e+05]); works after 2016b
set(gca,'XTick',[0:3e+05:27e+05]);
%xticklabels({'0', etc... --- necessary???

%name = strcat(sampleName,' (cutTot= ',num2str(cutTot),', cutRat= ',...
    %num2str(cutRat),')','-->',strcat(num2str(aboveCut),' insertions'));
%name = [name newline strcat(num2str(aboveCut),' insertions')]; for >=2017b
%title(name);
xlabel('chromosome coordinate');
ylabel('insertion sites / 10 kB');


%% for BEI_concat
%linear bar graph: coordinate vs. total reads, passing filters
figure;
edges = [1:10000:4920000];
hist(CTcut(:,1),edges); %,'FaceColor',[0.8 0.2 0,2]); 
xlim([-40000 4900000]);
%xticks([0:3e+05:27e+05]); works after 2016b
set(gca,'XTick',[0:3e+05:49e+05],'FontSize',13);
%xticklabels({'0', etc... --- necessary???

%name = strcat(sampleName,' (cutTot= ',num2str(cutTot),', cutRat= ',...
    %num2str(cutRat),')','-->',strcat(num2str(aboveCut),' insertions'));
%name = [name newline strcat(num2str(aboveCut),' insertions')]; for >=2017b
%title(name);
xlabel('chromosome coordinate');
ylabel('insertion sites / 10 kB');