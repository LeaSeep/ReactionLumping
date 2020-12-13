function [preproccessedModel]=ReadIn_iJ904()
    %%Load in model 
    fullModel=readCbModel('../data/iJR904.xml');
    
    %%
    [~,~,rxnsData]=xlsread('../data/ExampleExcelSheet_ECOLI.xlsx','Rxns');
    [~,~,metData]=xlsread('../data/ExampleExcelSheet_ECOLI.xlsx','Mets');
    rxnsData(1,:)=[];
    metData(1,:)=[];
    
    fullModel.mets2=replace(fullModel.mets,'[','_');
    fullModel.mets2=erase(fullModel.mets2,']');
    fullModel.mets2=replace(fullModel.mets2,'__','_');
    
    matchIDX=cellfun(@(s) find(strcmp(s,metData(:,1))),fullModel.mets2,'UniformOutput',false);
    matchIDX(find(cellfun('isempty',matchIDX)))={'Unknown'};
    matchIDX(find(cellfun(@(s) length(s)>1 ,matchIDX)))
    fullModel.DeltaG_m_std={};
    for i=1:length(matchIDX)
        if ~strcmp(metData(cell2mat(matchIDX(i)),2),'Unknown')
            fullModel.DeltaG_m_std_raw(i)=metData(cell2mat(matchIDX(i)),2);
            fullModel.pKa_offset(i)=metData(cell2mat(matchIDX(i)),8);
            fullModel.pH_offset(i)=metData(cell2mat(matchIDX(i)),9);
            fullModel.DeltaG_m_std(i)={cell2mat(fullModel.DeltaG_m_std_raw(i))+cell2mat(fullModel.pKa_offset(i))+cell2mat(fullModel.pH_offset(i))};
        else
            fullModel.DeltaG_m_std(i)={nan};
        end 
    end
    fullModel.DeltaG_m_std=fullModel.DeltaG_m_std;
    fullModel.DeltaG_m_std_match=fullModel.mets2;
    %% Manual curation
    % following also publication
    fullModel=removeRxns(fullModel,'CBLAT_DELETE');
    fullModel=removeRxns(fullModel,'CBIAT_DELETE');
    fullModel=removeRxns(fullModel,'BSORx');
    fullModel=removeRxns(fullModel,'BSORy');
    fullModel=removeRxns(fullModel,'DMSOR1');
    fullModel=removeRxns(fullModel,'DMSOR1e');
    fullModel=removeRxns(fullModel,'DMSOR2');
    fullModel=removeRxns(fullModel,'DMSOR2e');
    
    fullModel.rxns(find(strcmp(fullModel.rxns,'NACt')))={'NACUP'};
    fullModel.rxns(find(strcmp(fullModel.rxns,'ASPO3_1')))={'ASPO3'};
    fullModel.rxns(find(strcmp(fullModel.rxns,'ASPO4_1')))={'ASPO4'};
    fullModel.rxns(find(strcmp(fullModel.rxns,'ASPO5_1')))={'ASPO5'};
    fullModel.rxns(find(strcmp(fullModel.rxns,'ASPO6_1')))={'ASPO6'};
    
    % Remove copy entries always choose first copy
    fullModel.rxns(find(contains(fullModel.rxns,'_copy1')))=erase(fullModel.rxns(find(contains(fullModel.rxns,'_copy1'))),'_copy1');
    fullModel=removeRxns(fullModel,fullModel.rxns(find(contains(fullModel.rxns,'_copy'))));
    toRemove=[string('btnso[c]');string('dmso[c]');string('dmso[e]')];
    fullModel=removeMetabolites(fullModel,toRemove);
%% Output
     preproccessedModel=fullModel;
end