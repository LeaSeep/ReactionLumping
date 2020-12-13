function [alpha,L,remainAnIssue]=lumpReactions(model,idxToBeLumped,allowLumpingTransport)
%% The function lumpReaction takes a metabolic model as input as well as an indices-array
%  Indicating which metabolites have no available information for energies of formation. 
%  The output alpha matrix can be directly integrated in TMFA framework,whereas
%  the L output can be used to calculate the standard gibbs free energy fot
%  the new lumped reaction.
%% INPUT
%--model: metabolic model with at least the following fields:
%        -.S  
%        -.rxns & -.mets
%        -.lb & .-ub 
%--idxToBeLumped: indices of all metabolites that have no available
%                 information 
%--allowLumpingTransport: boolean, whether to allow transport reactions to
%                         be lumped or not
%% OUTPUT
%--alpha: matrix that has as many rowas as found lumped reactions and as
%many cols as there are (original) reactions, with a 1 indicating which
%original reactions are involved in the lumped reaction
%--L: is a model struc with the elements:
%       -.S     : as many cols as there Lumped reactions and lumped
%                 stoichiometric coefficients
%       -.mets  : all orignal metabolites
%       -.rxns  : names for lumped reactions
%--remainAnIssue: Metabolites that entered with unavailable energy
%                 information and were not eliminated with our procedure
%%
    disp("Lumping Procedure")
    %% Make sure to have splitted Model
    if(any(model.lb<0))
        splitModel=convertToIrreversible(model);
    else
         splitModel=model;
    end
    if(allowLumpingTransport==1)
        fprintf('Allow the Lumping of Transport Rxns: TRUE');
    end
    %% Identify reactions which are are presented to lumping procedure
    % Procedure to automatically exlucde biomass and export reactions and
    % transport reactions if corresponding paramtere was set
    % this can be skipped or expanded in order which reactions one wants to
    % be allowed to lump
    
    % For Mets: reduce all to without known Energy
    Aeq=splitModel.S(idxToBeLumped,:);
    [No.allMet,No.allrxns]=size(splitModel.S);
    % For Rxns: reduce all rxns to rxns involving Mets without unknwon E
    isInvolvingUnknown=zeros(size(Aeq,2),1);
    for i=1:size(Aeq,2)
        if(nnz(Aeq(:,i))>=1)
            isInvolvingUnknown(i)=1;
        end
    end
    % Reactions that are not allowed to lump are
    % BioMass, Exchange reactions and Transport reactions
    idxBio=find(contains(splitModel.rxns,'BIOMASS','IgnoreCase',true));
    idxEX=find(findExcRxns(splitModel));
    TransRxns=findTransRxns(splitModel);
    idxTrans=findRxnIDs(splitModel,TransRxns);
    if ismember(idxBio,find(isInvolvingUnknown)) 
        disp('Biomass is already a lumped reaction taking out !')
        isInvolvingUnknown(idxBio)=0;
    end
    if any(ismember(idxEX,find(isInvolvingUnknown)) )
        disp('Exchange Reactions are not to be lumped')
        isInvolvingUnknown(idxEX)=0;
    end
    if(~allowLumpingTransport)
        if any(ismember(idxTrans,find(isInvolvingUnknown)) )
            disp('Transport Reactions are not to be lumped - as setted Option')
            isInvolvingUnknown(idxTrans)=0;
        end
    end
    %
    fprintf('No. rxns involving met with unknown E before (splittedModel!): %d\n',nnz(isInvolvingUnknown))
    %% Reduce original Stochiometric prepped for Group Finding
    % reduce matrix to all meatbolites of unknwon E plus corresponding
    % reactions that are involving those
    G=Aeq(:,find(isInvolvingUnknown));
    row_all_zeros = find(all(G == 0,2));
    idxToBeLumped(row_all_zeros)=[];
    G=G(setdiff(1:size(G,1),row_all_zeros),:);
   
    %% Group Finding
    % find metabolites that participate in at least one reaction together
    % and indicate them as one group;(can of course consist only of a single metabolite)
    
    GroupOfMets=zeros(size(G,1),1);
    counter=0;
    arrayOfMets=[1:length(GroupOfMets)];
    while ~isempty(arrayOfMets)
        counter=counter+1;
        tmpRow=G(arrayOfMets(1),:);
        idxInvolvedRXNS=find(tmpRow~=0);
        idxInvoledMETS=find(sum(abs(G(:,idxInvolvedRXNS)),2)~=0);
        if(nnz(GroupOfMets(idxInvoledMETS))~=0)
            idxAlreadyDependent=find(GroupOfMets(idxInvoledMETS)~=0);
            tmp=GroupOfMets(idxInvoledMETS);
            for k=idxAlreadyDependent'
                dependentOn=tmp(k);
                idxToAdd=find(GroupOfMets==dependentOn);
                idxInvoledMETS=unique([idxInvoledMETS;idxToAdd]);
            end
        end
        GroupOfMets(idxInvoledMETS)=counter;
        arrayOfMets=setdiff(arrayOfMets,idxInvoledMETS);
    end
    %% Setting up variables
    nLumped=0;
    Lumped=[];%dim(No.allMet x No.Lumped)
    alpha=[]; %dim(No.Lumped x rxns)
    remainAnIssue=[];
    rxnsNames=[];
    [GC,GR] = groupcounts(GroupOfMets);
    [~,orderBasedOnGroupCount]=sort(GC,'descend');
    %% Get Reduced Matrix R
    % can be furhter reduced fitting personal needs
    R_r=model.S(:,find(isInvolvingUnknown)); % All metabolites (to get correct lumped reaction)
                                             % needed to fit dimensions as in here complete zero rows for
                                             % metabolites are still included
    reduceModel=removeRxns(model,model.rxns(~isInvolvingUnknown)); %all rxns involving met of unknownDeltaG + all present mets in those
    R=reduceModel.S;
    rel_idxToBeLumped=findMetIDs(reduceModel,model.mets(idxToBeLumped));
    fullModel_idx=findMetIDs(reduceModel,model.mets);                                         
    %% Lumping Procedure
    % implemented a feasibility program that aims to find a linear
    % combination of reactions that are involing mets with unknown E to
    % elimnate those
    timeCounter=1;
    textprogressbar('Lumping:')
    for i=GR(orderBasedOnGroupCount)'% sort after size of groups
        metsToLump_abs=idxToBeLumped(find(GroupOfMets==i)); % is absoluteIdx
        metsToLump_rel=findMetIDs(reduceModel,model.mets(metsToLump_abs));
        %integer Cut Matrix init - to avoid finding the same solution
        intCutMat=[];
        if (i>1) & ~isempty(alpha)
       	 intCutMat=[alpha(:,logical(isInvolvingUnknown))];
        end
        %% Group component
        [lumpVector_gr,lumpedRxn_gr,reactionsCoeff_gr]=basicLump(R,metsToLump_rel,intCutMat);
        timeCounter1=timeCounter1+1;
         
        if isempty(lumpVector_gr)
            Feasbility=0;
        else
           Feasbility=1;
        end
        if Feasbility==1
            %--Group Method was Feasible
            %get which rxns were used and put to alpha; dimension (lumpedRxns x rxns)
            %get the freshly lumped rxns and put to L; dimension (mets x No.Lumped)
            %lumpedRxn_gr 
            nLumped=nLumped+1;
            rxnsNames=[rxnsNames;cellstr(string(strcat('Lumped',string(nLumped))))];
            % Get Lumped and adjust to correct Dimensions
            tmp_Lumped=R_r*reactionsCoeff_gr;
            tmp_Lumped(abs(tmp_Lumped)<1e-5)=0;
            Lumped(:,nLumped)=tmp_Lumped;
            % Adjust alphaMat to correct Dimensions
            tmp_alphaRow=zeros(1,No.allrxns);
            absIDX=find(isInvolvingUnknown);
            relIDX=find(lumpVector_gr);
            LumpedRxnsIDX=absIDX(relIDX);
            tmp_alphaRow(LumpedRxnsIDX)=1;
            alpha(nLumped,:)=[tmp_alphaRow];
            %
            timeStamp=(timeCounter/length(GR))*100;
            textprogressbar(timeStamp)
        else
            %% Sequential Component is needed
            %Group Method could not find a solution for the group of mets
            %thats why we are SEQuentially checking the mets in the group
            if(size(metsToLump_abs,1)~=1)
                metsToLump_abs=metsToLump_abs';
            end
            tmp_timeCounter=1;
            fixedTimeStep=length(metsToLump_abs);
            while(~isempty(metsToLump_abs))
                timeStamp=((timeCounter-1)/length(GR))*100+(tmp_timeCounter/(fixedTimeStep*length(GR)))*100;
                textprogressbar(timeStamp)
                tmp_timeCounter=tmp_timeCounter+1;
                
                k=metsToLump_abs(1);
                AllMetToBeLumped_tmp=model.mets(metsToLump_abs);
                [~,noriginalRxns]=size(R);
                absMetIDX=k;
                relMetIDX=findMetIDs(reduceModel,model.mets(absMetIDX)); 
                
                %integer Cut Matrix init - to avoid finding the same solution
                intCutMat=[];
                if (i>1) & ~isempty(alpha)
                    intCutMat=[alpha(:,logical(isInvolvingUnknown))];
                end
                
                %Sequential Lumping program for each group metabolite
                [lumpVector_seq,lumpedRxn_seq,reactionsCoeff_seq]=basicLump(R,relMetIDX,intCutMat);
                             
                M=R; %M may be expanded with temporarly found Lumped rxns
                counter=1;
                LoopUpMatch={};
                LoopUpMatch_rxns=nan(size(R,1),1);
                
                if(~isempty(lumpVector_seq))
                    LoopUpMatch{counter}=[find(lumpVector_seq)];
                    LoopUpMatch_rxns(:,counter)=lumpedRxn_seq;
                else
                    remainAnIssue=[remainAnIssue,model.mets(absMetIDX)];
                    metsToLump_abs(1)=[];
                    continue
                end
                 
                while(length(find(lumpedRxn_seq(rel_idxToBeLumped)~=0))~=0  & ~isempty(lumpVector_seq))
                    % Sequential program from one group metabolite with
                    % expanding metabolites set that is constrained to 0
                    
                    leftOverMets= find(lumpedRxn_seq~=0); 
                    leftOverMets=intersect(leftOverMets,rel_idxToBeLumped);
                    relMetIDX=[relMetIDX leftOverMets(1)];
                    %expand M with latest solution
                    M(:,end+1)=lumpedRxn_seq;
                   
                    [lumpVector_seq,lumpedRxn_seq,reactionsCoeff_seq]=basicLump(M,relMetIDX,[]);
                    
                    if(isempty(lumpedRxn_seq)) %no solution can be found
                        break
                    end
                    if ~isempty(lumpedRxn_seq(rel_idxToBeLumped)) & nnz(lumpedRxn_seq(rel_idxToBeLumped))~=0
                        counter=counter+1;
                        LoopUpMatch{counter}=[find(lumpVector_seq)];
                        LoopUpMatch_rxns(:,counter)=lumpedRxn_seq; %multiplikative factor of rxns
                    end    
                end
                
                if(isempty(lumpedRxn_seq))
                    metsToLump_abs(1)=[];
                    remainAnIssue=[remainAnIssue,model.mets(absMetIDX)];
                    continue
                end
                if(sum(abs(lumpedRxn_seq(rel_idxToBeLumped)))~=0)
                    %Fail sequential lumping
                    metsToLump_abs(1)=[];
                    remainAnIssue=[remainAnIssue,model.mets(absMetIDX)];
                    
                else
                   %Success seuqntial Lumping
                   nLumped=nLumped+1;
                   rxnsNames=[rxnsNames;cellstr(string(strcat('Lumped',string(nLumped))))];
                   zeroMat=zeros(size(R_r,1),size(LoopUpMatch_rxns,2));
                   fullModel_idx(fullModel_idx==0)=[];
                   zeroMat(fullModel_idx,:)=LoopUpMatch_rxns;
                   LoopUpMatch_rxns=zeroMat;
                   
                   if sum(abs(LoopUpMatch_rxns(rel_idxToBeLumped,end)))==0
                        %first added reaction is already final lumped reaction
                        interMedLumped=0;
                   else
                        interMedLumped=LoopUpMatch_rxns;
                   end
                   
                   %Create Lumped reaction
                   if length(reactionsCoeff_seq)==size(R_r,2)
                       tmp_Lumped=[R_r]*reactionsCoeff_seq;
                       tmp_Lumped(abs(tmp_Lumped)<1e-5)=0;
                   else
                       tmp_Lumped=[R_r interMedLumped]*reactionsCoeff_seq;
                       tmp_Lumped(abs(tmp_Lumped)<1e-5)=0;
                   end
                   Lumped(:,nLumped)=tmp_Lumped;
                   tmp_Lumped(find(abs(tmp_Lumped)<1e-5))=0;
                   
                   %make Matrix alpha
                    if isempty(LoopUpMatch)
                        %addToAlpha adjust accordinlgy 
                        tmp_alphaRow=zeros(1,No.allrxns);
                        absIDX=find(isInvolvingUnknown);
                        relIDX=find(lumpVector_seq);
                        LumpedRxnsIDX=absIDX(relIDX);
                        tmp_alphaRow(LumpedRxnsIDX)=1;
                        alpha(nLumped,:)=[tmp_alphaRow];    
                    else
                        if sum(abs(lumpVector_seq(noriginalRxns:end)))~=0
                           lumpVector_seq(unique(vertcat(LoopUpMatch{:})))=1;
                        end
                        lumpVector_seq(noriginalRxns+1:end)=[];
                        %addToAlpha adjust accordinlgy 
                        tmp_alphaRow=zeros(1,No.allrxns);
                        absIDX=find(isInvolvingUnknown);
                        relIDX=find(lumpVector_seq);
                        LumpedRxnsIDX=absIDX(relIDX);
                        tmp_alphaRow(LumpedRxnsIDX)=1;
                        alpha(nLumped,:)=[tmp_alphaRow]; 
                    end    
                end
                tmp=setdiff(AllMetToBeLumped_tmp,reduceModel.mets(relMetIDX));
                metsToLump_abs=findMetIDs(model,tmp);  
            end
        end
    end
    %% Summary Numbers
    textprogressbar(100);
    textprogressbar('done');
    if ~isempty(remainAnIssue)
        lumpedMets=length(setdiff(model.mets(idxToBeLumped),unique(remainAnIssue)));

        fprintf('Resulting in %d Lumped Reactions\nEqualing out %d mets with unknown Energy.\n',nLumped,lumpedMets);
    end
    %% Remove duplicate rows
    before=size(alpha,1);
    [~, alpha_tmp] = uniquetol(alpha, 'byrows', true);
    alpha=alpha(sort(alpha_tmp), :);      %the unique rows, in the original order
    after=size(alpha_tmp,1);
    doubleEntries=before-after;
    fprintf('Removed %d doubled entries\n',doubleEntries)
    %% Statisics
    allZerosCols = find(sum(abs(alpha(:,logical(isInvolvingUnknown)))) == 0);
    %% Output
    if(nLumped==0)
        alpha=[];
        L.S=[];
        L.mets=model.mets;
        L.rxns=[];
    else
        alpha=alpha;
        L.S=Lumped;
        L.mets=model.mets;
        L.rxns=rxnsNames;
    end
end
