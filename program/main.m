%% SetUp the workspace
% Load CobraToolbox; code was test with version: 2.29.2
initCobraToolbox()
%% Read In Data
% Needed is Metabolic-Network model and table with available thermodynamic
% information (Gibbs free Energy of Formation), this may be substitute with
% a indices of metabolites with unavailable information

% Read in Model

%Ecoli 
fullModel=ReadIn_iJ904();
%Bacillus subitils
%fullModel=ReadIn_Bacillus();
%fullModel.description='Bacillus';
%redHuman
%fullModel=ReadIn_redHuman();
%fullModel.description='redHuman';
% Make an irreversible Model
fullModel=convertToIrreversible(fullModel);
% Identify metabolite with missing change of Gibbs free energy of formation
% Specifiy metabolites with missing data a indices list of metabolites with
% missing Gibbs free Energy of Formation is needed.
% If available thermodynamic information is given, like in our
% examples the following retrieves the indices 
% (please see the corresponding ReadIn_model function to see the inclusion 
% of excel sheets of thermodynamic information and corresponding models)
if iscell(fullModel.DeltaG_m_std)
   fullModel.DeltaG_m_std=cell2mat(fullModel.DeltaG_m_std);
end
idx.toBelumped=find((isnan(fullModel.DeltaG_m_std)==1));

%% Lumping procedure
[alpha,lumpedModel,metsToDelete]=lumpReactions(fullModel,idx.toBelumped,true);

%% saving
match='.';
FileName=str2mat(eraseBetween(fullModel.description,match,fullModel.description(end),'Boundaries','inclusive'));
FileName1=['../result/' FileName '_LumpedModel.mat'];
save(FileName1, '-struct', 'lumpedModel');
FileName2=['../result/' FileName '_alphaMatrix.mat'];
save(FileName2, 'alpha');
%% DONE
%% This section is to retrieve the timing data for the paper figures 2-3
% Lumping with timings of different components of the procedure
[alpha ,lumpedModel ,metsToDelete,Timing]=lumpReactions_timing(fullModel,idx.toBelumped,true);

%% saveData for Figure 2
dataToSave={'before' 'known'   'rxns' (Timing.allRxns-Timing.NoRxnsUnconstrained_before_all)/Timing.allRxns;
            'before' 'unknown' 'rxns' (Timing.NoRxnsUnconstrained_before_all)/Timing.allRxns                   ;
            'after'  'known'   'rxns' (Timing.allRxns-Timing.NoRxnsUnconstrained_after)/Timing.allRxns     ;
            'after'  'unknown' 'rxns' (Timing.NoRxnsUnconstrained_after)/Timing.allRxns     ;
            'before' 'unknown' 'mets' Timing.NoMetsUnknown/Timing.allMets;
            'after' 'success' 'mets'  Timing.NoMetsCouldBeEliminated_after/Timing.allMets;
            'after' 'remain' 'mets' (Timing.NoMetsUnknown-Timing.NoMetsCouldBeEliminated_after)/Timing.allMets};
% Convert cell to a table 
T = cell2table(dataToSave(1:end,:));
FileName=str2mat(eraseBetween(fullModel.description,match,fullModel.description(end),'Boundaries','inclusive'));
FileName=['../result/' FileName 'Figure1Data.csv'];
% Write the table to a CSV file
writetable(T,FileName)

%% saveData for Figure 3
%Data_part1
toTable.time=[Timing.group(:,2);Timing.initialLumping';Timing.nestedWhileLoop(:,3)];
toTable.type=[string(repelem('g',size(Timing.group,1))');
            string(repelem('s',size(Timing.initialLumping,2))');
            string(repelem('s',size(Timing.nestedWhileLoop,1))')];
T = struct2table(toTable);

FileName=str2mat(eraseBetween(fullModel.description,match,fullModel.description(end),'Boundaries','inclusive'));
FileName=['../result/' FileName 'Figure2Data_part1.csv'];
% Write the table to a CSV file
writetable(T,FileName)

%Data_part2
dataToSaveGroup_NoMembers.T=[Timing.Group_NoMembers];
T = struct2table(dataToSaveGroup_NoMembers);
FileName=str2mat(eraseBetween(fullModel.description,match,fullModel.description(end),'Boundaries','inclusive'));
FileName=['../result/' FileName 'Figure2Data_part2.csv'];
% Write the table to a CSV file
writetable(T,FileName)

%Data3
dataToSave_MetsLeft.MetsLeft=[Timing.MetsLeft_reason{:,1}]';
dataToSave_MetsLeft.Type={Timing.MetsLeft_reason{:,2}}';
T = struct2table(dataToSave_MetsLeft);
FileName=str2mat(eraseBetween(fullModel.description,match,fullModel.description(end),'Boundaries','inclusive'));
FileName=['../result/' FileName 'Figure2Data_part3.csv'];
% Write the table to a CSV file
writetable(T,FileName)