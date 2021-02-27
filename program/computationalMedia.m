function outputModel=computationalMedia(model,method)
%% Options are
% Reed Table S6
%Feist et al , s2 in silico (computational minmal media

if contains(method,'Reed_aerobicGlucose')
    %taken from Table S6 Reed
    disp('Using Reed minimal Media -aerobic glucose Growth - ECOLI')
    %Fixed
    model.lb(contains(model.rxns,'EX_'))=0;
    model.lb(find(contains(model.rxns,'EX_nh4')))=-1000;%
    model.lb(find(contains(model.rxns,'EX_so4')))=-1000;%
    model.lb(find(contains(model.rxns,'EX_pi')))=-1000;%
    model.lb(find(contains(model.rxns,'EX_co2')))=-1000;%
    model.lb(find(contains(model.rxns,'EX_h_')))=-1000;%
    model.lb(find(contains(model.rxns,'EX_h2o')))=-1000;%
    model.lb(find(contains(model.rxns,'EX_fe2')))=-1000;% 
    model.lb(find(contains(model.rxns,'EX_k')))=-1000; %
    model.lb(find(contains(model.rxns,'EX_na1')))=-1000;%
    % vary depedning on simulation
    model.lb(find(contains(model.rxns,'EX_o2')))=-1000;
    if(~isempty(find(contains(model.rxns,'EX_glc__D_e'))))
         model.lb(find(contains(model.rxns,'EX_glc__D_e')))=-10;
    else
         model.lb(find(contains(model.rxns,'EX_glc__D_e_r')))=-10;
    end
    model.lb(find(contains(model.rxns,'EX_glyc_')))=0;
elseif contains(method,'Henry_Glucose')
     %taken from Table S6 Reed
    disp('Using Henry Media Formulation of Glucose Growth - BACILLUS')
    model.lb(contains(model.rxns,'EX_'))=0;
    model.lb(find(contains(model.rxns,'EX_cpd00027_e')))=-1000;
    model.lb(find(contains(model.rxns,'EX_cpd00013_e')))=-1000;
    model.lb(find(contains(model.rxns,'EX_cpd00009_e')))=-1000;
    model.lb(find(contains(model.rxns,'EX_cpd00048_e')))=-1000;
    model.lb(find(contains(model.rxns,'EX_cpd00063_e')))=-1000;
    model.lb(find(contains(model.rxns,'EX_cpd00011_e')))=-1000;
    model.lb(find(contains(model.rxns,'EX_cpd10516_e')))=-1000; 
    model.lb(find(contains(model.rxns,'EX_cpd00067_e')))=-1000;
    model.lb(find(contains(model.rxns,'EX_cpd00001_e')))=-1000;
    model.lb(find(contains(model.rxns,'EX_cpd00205_e')))=-1000;
    model.lb(find(contains(model.rxns,'EX_cpd00254_e')))=-1000;
    model.lb(find(contains(model.rxns,'EX_cpd00971_e')))=-1000;
    model.lb(find(contains(model.rxns,'EX_cpd00007_e')))=-1000;
    model.lb(find(contains(model.rxns,'EX_cpd00099_e')))=-1000;
    model.lb(find(contains(model.rxns,'EX_cpd00058_e')))=-1000;
    model.lb(find(contains(model.rxns,'EX_cpd00149_e')))=-1000;
    model.lb(find(contains(model.rxns,'EX_cpd00030_e')))=-1000;
    model.lb(find(contains(model.rxns,'EX_cpd00034_e')))=-1000;  
end
 
outputModel=model;
    
end