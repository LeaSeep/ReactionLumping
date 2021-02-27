%% Calculate the \DeltaG standard for a reaction
%% Documentation
% INPUT:
%     model (irreversible!)
%           classical SBML model; additional fields
%           .DeltaG_m_std
%               for each Metabolite it deltaG_formation^0 value, missing
%               values should be nan 
% OUTPUT:
%     resultModel
%        with a new field 'DeltaG_r_total'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [resultModel]=calcDeltaG_r(model)
%% Check if reaction involves species from different compartments
% those need an additional adjustment 
 DeltaG_rxns_std=nan(length(model.rxns),1);
 isCellularOnly=logical(~findExcRxns(model));
 if(strcmp('cell',class(model.DeltaG_m_std)))
	 model.DeltaG_m_std=cell2mat(model.DeltaG_m_std);
 end
%% Filter out all reactions involving mets with unknown DeltaG_m
 RxnsWithUnknownMets=findRxnsFromMets(model,model.mets(find(isnan(model.DeltaG_m_std))));
 fprintf('%d reaction(s) which Delta G_r_total cannot be calculated\n',length(RxnsWithUnknownMets))
 DeltaG_rxns_std(findRxnIDs(model,RxnsWithUnknownMets))=NaN;
 %reduce stoich and model.DeltaG_m_std as matlab does not deal well with
 %Nan included in the DeltaG_m_std vector -> would result in all Nan
 tmp=model;
 if ~isfield(model,'pH_offset') 
     model.pH_offset={repmat(0,length(model.mets),1)};
 else
     % Ecoli data specific as not already incorporated in given DeltaGm values
     model.pH_offset(cellfun(@isempty, model.pH_offset))={0};
 end
 DeltaG_stt_intracellular=(model.DeltaG_m_std)'-cell2mat(model.pH_offset)';
 
 %change all NaN within DeltaG_stt_intracellular to 0 ; do the calculcation
 %change all DeltaG_rxns involving mets of unknown to NaN
 DeltaG_stt_intracellular(isnan(DeltaG_stt_intracellular))=0;
 if(size(DeltaG_stt_intracellular,1)<size(DeltaG_stt_intracellular,2))
    DeltaG_stt_intracellular=DeltaG_stt_intracellular';
 end
 DeltaG_rxns_intracellular=tmp.S'*DeltaG_stt_intracellular;
 DeltaG_rxns_intracellular(findRxnIDs(model,findRxnsFromMets(model,model.mets(isnan(model.DeltaG_m_std)))))=nan;
 
 DeltaG_rxns_transport=zeros(length(DeltaG_rxns_intracellular),1);
 %%if reaction is transport reaction energy term of transport must be added
 for idxTransport=find(~isCellularOnly)'
     all_metsInvolved=findMetsFromRxns(model,model.rxns(idxTransport));
     exp_mets=all_metsInvolved(model.compartment(findMetIDs(model,all_metsInvolved))=='e');
     DeltaTransport=0;
     for everyExp=1:length(exp_mets)
         %erase CompInfo to only have Met
         tmpMet_e=exp_mets(everyExp);
         tmpMet_raw=eraseBetween(tmpMet_e,'[',']','Boundaries','inclusive');
          
         netCharge=0;
         transportH=0;
         if(strcmp('h',tmpMet_raw))
             transportH=abs(model.S(findMetIDs(model,tmpMet_e),idxTransport));
         else
           netCharge=model.metCharges(findMetIDs(model,tmpMet_e))*abs(model.S(findMetIDs(model,tmpMet_e),idxTransport));  
         end
         DeltaTransport(everyExp)=netCharge*constant.deltaPsi*constant.F-2.3*transportH*constant.R_kcal*constant.T0;
     end
     DeltaGTransport_total=sum(DeltaTransport);
     DeltaG_rxns_transport(idxTransport)=DeltaGTransport_total;
 end
 
 DeltaG_rxns_intracellular(findRxnIDs(model,RxnsWithUnknownMets))=NaN;
 
 %% Output
 resultModel=model;
 resultModel.DeltaG_rxns_intracellular=DeltaG_rxns_intracellular;
 resultModel.DeltaG_rxns_transport=DeltaG_rxns_transport;
 resultModel.DeltaG_r_total=DeltaG_rxns_intracellular+DeltaG_rxns_transport;
 
end