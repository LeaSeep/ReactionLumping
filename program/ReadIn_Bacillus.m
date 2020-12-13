function [model]=ReadIn_Bacillus()
model=readCbModel('../data/Bacillus.xml');
disp('Read in Bacillus')
model.mets=eraseBetween(model.mets,'[',']','Boundaries','inclusive');
metData=load("../data/BacillusMetData.mat");
metData=metData.BacillusDataMediumSIzeNetwork;
rxnsData=load("../data/BacillusRxnsData.mat");
rxnsData=rxnsData.BacillusDataMediumSIzeNetworkS1;
rxnsData(contains(rxnsData(:,1),'BIOMASS'),:)=[];

rxnsData(1,:)=[];  %'REACTION ID' 'DELTAG (kcal/mol)'
metData(1,:)
metData(1,:)=[];   %'COMPOUND ID' 'DELTAG (kcal/mol)'

metData(cellfun(@(s) strcmp(s,'UNKNOWN'), metData(:,3)),:)=[];

orderRxnsData=findRxnIDs(model,rxnsData(:,1)); % there is one idx= 0 (belonging to Biomass Rxns)
orderMetsData=findMetIDs(model,metData(:,1));
metData(find(orderMetsData==0),:)=[];
metData(strcmp(metData(:,1),'cpd00464'),:)=[];
orderMetsData=findMetIDs(model,metData(:,1));

model.DeltaG_m_std=nan(length(model.mets),1);
model.DeltaG_m_std(orderMetsData)=cell2mat(metData(:,3));

fields  = fieldnames(model);
mask=zeros(length(fields),1);
mask(find(strcmp(fields, 'S')))=1;
mask(find(strcmp(fields, 'rxns')))=1;
mask(find(strcmp(fields, 'mets')))=1;
mask(find(strcmp(fields, 'DeltaG_m_std')))=1;
mask(find(strcmp(fields, 'lb')))=1;
mask(find(strcmp(fields, 'ub')))=1;
mask(find(strcmp(fields, 'c')))=1;
model = rmfield(model,fields(~mask));


end