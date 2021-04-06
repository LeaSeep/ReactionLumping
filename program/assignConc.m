function [result]=assignConc(model,type)
 disp('ASIGN CONCENTRATION')
%This numbers are taken from  Table S6
%Hamilton JJ, Dwivedi V, Reed JL.
%Quantitative assessment of thermodynamic constraints on the solution space of genome-scale metabolic models

%% Get compartment data from .mets 
model.metCompSymbol=regexprep(model.mets, '^.*[', '', 'lineanchors');
model.metCompSymbol=regexprep(model.metCompSymbol, ']', '', 'lineanchors');
model.old_mets=model.mets;
model.mets=eraseBetween(model.mets,'[',']','Boundaries','inclusive');
%% General boundaries for internal and external Concentration

conc_ranges.lb=log(ones(length(model.mets),1)*10^(-3)*10^-3); %M  unit
conc_ranges.ub=log(ones(length(model.mets),1)*20*10^-3);      %M  unit

%% Table
%extra cellular
if type=="ECOLI"
if(true)
    conc_ranges.lb(find((strcmp(model.mets,'pi') +strcmp(model.metCompSymbol,'e'))==2))=log(2*10^-3);
    conc_ranges.ub(find((strcmp(model.mets,'pi') +strcmp(model.metCompSymbol,'e'))==2))=log(2*10^-3);
    conc_ranges.lb(find((strcmp(model.mets,'so4') +strcmp(model.metCompSymbol,'e'))==2))=log(0.286*10^-3);
    conc_ranges.ub(find((strcmp(model.mets,'so4') +strcmp(model.metCompSymbol,'e'))==2))=log(0.286*10^-3);
    conc_ranges.lb(find((strcmp(model.mets,'nh4') +strcmp(model.metCompSymbol,'e'))==2))=log(9.52*10^-3);
    conc_ranges.ub(find((strcmp(model.mets,'nh4') +strcmp(model.metCompSymbol,'e'))==2))=log(9.52*10^-3);
    conc_ranges.lb(find((strcmp(model.mets,'na1') +strcmp(model.metCompSymbol,'e'))==2))=log(50*10^-3);
    conc_ranges.ub(find((strcmp(model.mets,'na1') +strcmp(model.metCompSymbol,'e'))==2))=log(50*10^-3);
    conc_ranges.lb(find((strcmp(model.mets,'k') +strcmp(model.metCompSymbol,'e'))==2))=log(3.19*10^-3);
    conc_ranges.ub(find((strcmp(model.mets,'k') +strcmp(model.metCompSymbol,'e'))==2))=log(3.19*10^-3);
    conc_ranges.lb(find((strcmp(model.mets,'fe2') +strcmp(model.metCompSymbol,'e'))==2))=log(0.01*10^-3);
    conc_ranges.ub(find((strcmp(model.mets,'fe2') +strcmp(model.metCompSymbol,'e'))==2))=log(0.01*10^-3);
    conc_ranges.lb(find((strcmp(model.mets,'co2') +strcmp(model.metCompSymbol,'e'))==2))=log(0.01*10^-3);
    conc_ranges.ub(find((strcmp(model.mets,'co2') +strcmp(model.metCompSymbol,'e'))==2))=log(0.01*10^-3);
    conc_ranges.lb(find((strcmp(model.mets,'o2') +strcmp(model.metCompSymbol,'e'))==2))=log(0.0082*10^-3);
    conc_ranges.ub(find((strcmp(model.mets,'o2') +strcmp(model.metCompSymbol,'e'))==2))=log(0.0082*10^-3);
    conc_ranges.lb(find((strcmp(model.mets,'glyc') +strcmp(model.metCompSymbol,'e'))==2))=0;
    conc_ranges.ub(find((strcmp(model.mets,'glyc') +strcmp(model.metCompSymbol,'e'))==2))=0;
    conc_ranges.lb(find((strcmp(model.mets,'glc__D') +strcmp(model.metCompSymbol,'e'))==2))=log(0.2*10^-3);
    conc_ranges.ub(find((strcmp(model.mets,'glc__D') +strcmp(model.metCompSymbol,'e'))==2))=log(0.2*10^-3);
    conc_ranges.lb(find((strcmp(model.mets,'h') +strcmp(model.metCompSymbol,'e'))==2))=log(0.00004*10^-3); %H+ extracellular
    conc_ranges.ub(find((strcmp(model.mets,'h') +strcmp(model.metCompSymbol,'e'))==2))=log(0.00004*10^-3); %H+ extracellular
end

if(true)
%intra cellular
    conc_ranges.lb(find((strcmp(model.mets,'h') + ~strcmp(model.metCompSymbol,'e'))==2))=log(0.0001*10^-3);
    conc_ranges.ub(find((strcmp(model.mets,'h') + ~strcmp(model.metCompSymbol,'e'))==2))=log(0.0001*10^-3);
    conc_ranges.ub(find((strcmp(model.mets,'co2') + ~strcmp(model.metCompSymbol,'e'))==2))=log(0.01*10^-3);
    conc_ranges.ub(find((strcmp(model.mets,'o2') + ~strcmp(model.metCompSymbol,'e'))==2))=log(0.0082*10^-3); 
    %conc_ranges.lb(find((strcmp(model.mets,'o2') + ~strcmp(model.metCompSymbol,'e'))==2))=log(10^-12); 
    conc_ranges.lb(find((strcmp(model.mets,'atp') + ~strcmp(model.metCompSymbol,'e'))==2))=log(8.13*10^-3);
    conc_ranges.ub(find((strcmp(model.mets,'atp') + ~strcmp(model.metCompSymbol,'e'))==2))=log(11.4*10^-3);
    conc_ranges.lb(find((strcmp(model.mets,'adp') + ~strcmp(model.metCompSymbol,'e'))==2))=log(0.437*10^-3); 
    conc_ranges.ub(find((strcmp(model.mets,'adp') + ~strcmp(model.metCompSymbol,'e'))==2))=log(0.704*10^-3);
end
if(false) %Maranas ranges for gases
    conc_ranges.lb(find((strcmp(model.mets,'h') + ~strcmp(model.metCompSymbol,'e'))==2))=log(10^-8);
    conc_ranges.ub(find((strcmp(model.mets,'h') + ~strcmp(model.metCompSymbol,'e'))==2))=log(0.000034);
    conc_ranges.ub(find((strcmp(model.mets,'co2') + ~strcmp(model.metCompSymbol,'e'))==2))=log(0.0014); 
    conc_ranges.lb(find((strcmp(model.mets,'co2') + ~strcmp(model.metCompSymbol,'e'))==2))=log(10^-8); 
    conc_ranges.ub(find((strcmp(model.mets,'o2') + ~strcmp(model.metCompSymbol,'e'))==2))=log(0.000055); 
    conc_ranges.lb(find((strcmp(model.mets,'o2') + ~strcmp(model.metCompSymbol,'e'))==2))=log(10^-8); 

    conc_ranges.lb(find((strcmp(model.mets,'atp') + ~strcmp(model.metCompSymbol,'e'))==2))=log(8.13*10^-3);
    conc_ranges.ub(find((strcmp(model.mets,'atp') + ~strcmp(model.metCompSymbol,'e'))==2))=log(11.4*10^-3); 
    conc_ranges.lb(find((strcmp(model.mets,'adp') + ~strcmp(model.metCompSymbol,'e'))==2))=log(0.437*10^-3);
    conc_ranges.ub(find((strcmp(model.mets,'adp') + ~strcmp(model.metCompSymbol,'e'))==2))=log(0.704*10^-3); 
     
end
end
%% Output
%change back to original mets (again with compartment info)
model.mets=model.old_mets;

result=model;
result.compartment=cell2mat(result.metCompSymbol);
result.conc_ub=conc_ranges.ub;
result.conc_lb=conc_ranges.lb;

end