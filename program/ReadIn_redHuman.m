function [model]=ReadIn_redHuman()
load('D:\Uni\Potsdam\4.Semester\ProjectWork\ProjectWork\data\additionalData\recon3_redHUMAN_curated.mat')
load('..\data\recon3_redHUMAN_curated.mat')
disp('Read in recon3 RedHuman model')
model.DeltaG_m_std=model.metDeltaGFstd;

model.DeltaG_m_std(model.DeltaG_m_std==10000000)=NaN;
fields  = fieldnames(model);
% Note: masking not necassary but cleaner as there are great number of unused
% fileds otherwise
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