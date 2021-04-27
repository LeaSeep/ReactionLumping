%% Perform Variability analysis of the DeltaGs before and after lumping
% choose model: Ecoli or Bacillus
clear
modelChoice="Ecoli";
%% SetUp the workspace
% Load CobraToolbox; code was test with version: 2.29.2
% initCobraToolbox(false)
cd("..\program")
%% Read in Model
if(modelChoice=="Ecoli")
    fprintf("Choosen model: 'Ecoli'\n")
    fullModel=ReadIn_iJ904();
    fullModel=computationalMedia(fullModel,'Reed_aerobicGlucose');  
    bigM=10^9;
elseif modelChoice=="Bacillus"
    fprintf("Choosen model: 'Bacillus'\n")
    fullModel=ReadIn_Bacillus(); 
    fullModel.description='Bacillus';
    fullModel=computationalMedia(fullModel,'Henry_Glucose');
    bigM=10^12;
else 
    error("Choose either 'Ecoli' or 'Bacillus'")
end
% get optimal value based on preassumed directionalities and specified
% computational Media
optVal=optimizeCbModel(fullModel);
optVal=optVal.f;
% choose level of optimality
percentage =0.9;
%% Model Settings
%Set boundaries to -1000 and 1000
%Ensure that the model has no preassumed directionalities for intracellular Rxns
%Allow all exchange reactions to excreed metabolites to the environment
%Reomve any blocked reactions left
fullModel.lb(fullModel.lb<-1000)=-1000;
fullModel.ub(fullModel.ub>1000)=1000;
fullModel.lb(findRxnIDs(fullModel,findIntracellularRxns(fullModel)))=-1000;
fullModel.ub(findRxnIDs(fullModel,findIntracellularRxns(fullModel)))=1000;
fullModel.ub(findExcRxns(fullModel)==1)=1000;
fullModel.lb(fullModel.c==1)=0; %BioMass
fullModel.ub(findRxnIDs(fullModel,{'ATPM'}))=7.6;%ATPM
fullModel.lb(findRxnIDs(fullModel,{'ATPM'}))=7.6;%ATPM
%Assign all reactions to be reversible
blockedRxns=findBlockedReaction(fullModel);
fullModel.lb(findRxnIDs(fullModel,blockedRxns))=0;
fullModel.ub(findRxnIDs(fullModel,blockedRxns))=0;
original=fullModel;

%% Assign Concentration for Extra and intracellular Concentration
if(modelChoice=="Ecoli")
    fullModel=computationalMedia(fullModel,'Reed_aerobicGlucose');
    fullModel=assignConc(fullModel,'ECOLI');
elseif modelChoice=="Bacillus"
    fullModel=computationalMedia(fullModel,'Henry_Glucose');
    fullModel=assignConc(fullModel,'Bacillus'); % Note: all metabolites get default ranges advised
end
%% Make an irreversible Model
fullModel=convertToIrreversible(fullModel,'orderReactions',true);

%% Calculate DeltaG_std for reactions based on given DeltaG_m
fullModel=calcDeltaG_r(fullModel); 
%% Identify metabolites without data -> propose to lumping
idx.toBelumped=find((isnan(fullModel.DeltaG_m_std)==1));
%% Run TMFA based model to get Concentration Boundaries
% The program getBoundaries is TMFA based while allowing violation of given
% concentration boundatiers, which is subjected to minimizatoin of such
% violation. All done while assuring TMFA constraints as well as an
% optimality constraint (based on FBA value)
% This is done for the model with and without Lumping
[conc_boundaries_withLump,x0_withLump]=getBoundaries(fullModel,bigM,idx.toBelumped,optVal,percentage);
%[conc_boundaries_withoutLump,x0_withoutLump]=getBoundaries(fullModel,bigM,[],optVal,percentage);

%% Run TMFA & Variability analysis for Delta G
for type=["withLumping","withOutLumping"]
    fprintf('Variability Analysis is running for\n')
    type
    if all(type=='withLumping')
        fullModel.conc_lb=conc_boundaries_withLump.lb;
        fullModel.conc_ub=conc_boundaries_withLump.ub;
        [Distr, optValTMFA, STAT, usedModel]=TMFA_LS(fullModel,bigM,idx.toBelumped,x0_withLump,optVal,percentage);
        rangesTable=myDeltaG_VA(usedModel,optValTMFA);
        savingName1=cell2mat({'../result/' str2mat(modelChoice) '_' str2mat(type) '_' 'Ranges_alternativeSolution.csv'});
        savingName2=cell2mat({'../result/' str2mat(modelChoice) '_' str2mat(type) '_' 'RxnsNames_alternativeSolution.csv'});
        writematrix(rangesTable,savingName1)
        writecell(fullModel.rxns,savingName2)
    else
        fullModel.conc_lb=conc_boundaries_withoutLump.lb;
        fullModel.conc_ub=conc_boundaries_withoutLump.ub;
        [Distr, optValTMFA, STAT, usedModel]=TMFA_LS(fullModel,bigM,[],x0_withoutLump,optVal,percentage);
        rangesTable=myDeltaG_VA(usedModel,optValTMFA);
        savingName1=cell2mat({'../result/' str2mat(modelChoice) '_' str2mat(type) '_' 'Ranges_alternativeSolution.csv'});
        savingName2=cell2mat({'../result/' str2mat(modelChoice) '_' str2mat(type) '_' 'RxnsNames_alternativeSolution.csv'});
        writematrix(rangesTable,savingName1)
        writecell(fullModel.rxns,savingName2)
    end
end
