%% Investigate into soft boundaries for metabolite concentrations
%% Documentation
% INPUT:
%     model (irreversible!)
%           classical SBML model; additional fields
%           .DeltaG_r_total : 
%              for each Reaction it deltaG^0 value, missing values should be nan
%           .DeltaG_m_std
%               for each Metabolite it deltaG_formation^0 value, missing
%               values should be nan (needed only if lumping should be
%               performed)
%           .conc_lb
%           .conc_ub
%               for each metabolite its boundaries (Literature says in
%               general for internal mets 10^-3:20mM 
%     bigM
%           high constant to be used
%     idxtoBelumped 
%           idx of metabolites which DeltaGm is missing (can be provided
%           empty to skip lumping
%     optVal: optimal Value of objective (eg Biomass) determined with FBA
%     percentage: level which is to be at least ensured (expected value of
%                 eg growth = percentage*optVal)
% OUTPUT:
%     metsConcRanges
%        determined concentration ranges
%     x0
%        intial solution to be given to TMFA as startpoint
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[metsConcRanges,x0]=getBoundaries(prepped_model,bigM,idxtoBelumped,optVal,percentage)
%% Manual add mets which are not allowed to violate their given boundaries
% Descicions made based on biological considerations (eg thight regulation
% of Protons within the cell)
notToLoose={'h[c]','pi[e]','so4[e]','nh4[e]','na1[e]','k[e]','fe2[e]','co2[e]','o2[e]','glyc[e]','glc__D[e]','atp[c]','adp[c]'};
%notToLoose={'g3pe[c]','g3pg[c]','h2o[c]','cmp[c]','trdox[c]','nadh[c]'};
%idx.essentAssocMets=[findMetIDs(prepped_model,notToLoose)];
notToLooseIDX=[findMetIDs(prepped_model,notToLoose)];
idx.essentAssocMets=setdiff(1:756,sort(notToLooseIDX));

n.essen=length(idx.essentAssocMets);
%% Parameter Setting 
original=prepped_model;
[n.Species,n.Rxns]=size(prepped_model.S);
R=constant.R_kcal;          
T=constant.T0;                 %[K] Temperature given
K=bigM;                        %upper boundary on DeltaG
deltaPsi=constant.deltaPsi;    %[mV]
deltapH=constant.deltapH;

%% Lumping
% search for metabolites with unknown Energy as they cause reactions to
% have undetermined deltaG values
% if idxtoBelumped is empty ('[]'), no lumping will be considered
  if ~isempty(idxtoBelumped)
        if iscell(prepped_model.DeltaG_m_std)
            prepped_model.DeltaG_m_std=cell2mat(prepped_model.DeltaG_m_std);
        end
        
        [alpha,lumpedModel,~]=lumpReactions(prepped_model,idxtoBelumped,1);
        
        % Update lumped Model with metNames + Energy
        lumpedModel.metNames=prepped_model.metNames;
        lumpedModel.DeltaG_m_std=prepped_model.DeltaG_m_std;
        if(isfield(prepped_model,'pH_offset')) %needed?
            lumpedModel.pH_offset=prepped_model.pH_offset;
        end
        if(isfield(prepped_model,'metCharges')) %needed?
            lumpedModel.metCharges=prepped_model.metCharges;
        end
        n.lumped=size(alpha,1);
        a=alpha;
        L=lumpedModel.S;
        n.lumped=size(alpha,1); 
        
        lumpedModel.lb=zeros(n.lumped,1);
        lumpedModel.ub=ones(n.lumped,1)*max(prepped_model.ub);
    else
        idx.noDeltaE_Rxns_noLumped=[];
        L=[];
        a=[];
        lumpedModel.DeltaG_r_std=[];
  end
       
%% Assinging deltaG to lumped reactions
    if ~isempty(L)
        disp('Assinging deltaG to Lumped Rxns')
        rxnIDXtoAssign=1:size(lumpedModel.S,2);
        lumpedModel.unKnownG_m=isnan(lumpedModel.DeltaG_m_std);
        %init result vector all NaN (for lumping reactions no NaN should be left after running calcDeltaG_r!)
        lumpedModel.DeltaG_r_std=nan(n.lumped,1);
        lumpedModel=calcDeltaG_r(lumpedModel);   
    else
         lumpedModel.DeltaG_r_total=[];
    end
%% Find out reactions that are not constraint due to missing deltaG_r and not included in any lumped reactions 
% this reactions are allowed to vary freely with respect to their deltaG
    idx.noDeltaE_Rxns=find(isnan(prepped_model.DeltaG_r_total));
    idx.noDeltaE_Rxns=[idx.noDeltaE_Rxns; find(prepped_model.c==1)];
    if ~isempty(a)
        idx.Rxns_inLumped=find(sum(a,1)==1)'; % included in any lumped reaction
    end
    if ~isfield(idx,'Rxns_inLumped') 
       idx.noDeltaE_Rxns_noLumped=idx.noDeltaE_Rxns;
    else
       idx.noDeltaE_Rxns_noLumped=setdiff(idx.noDeltaE_Rxns,idx.Rxns_inLumped);
    end
    n.noDeltaE_Rxns_noLumped=length(idx.noDeltaE_Rxns_noLumped);

%% Preparation of TMFA Matrices
%specified order of x
%x= [  v  deltaGr  ln(conc)  binary(lumped) binary_z(all) ]
    
    S=prepped_model.S;
    [n.Species,n.Rxns]=size(S);
    obj=prepped_model.c;
    disp('Objective Function to be maximized is:');
    idx.Bio=find(prepped_model.c==1);
    disp(prepped_model.rxns(idx.Bio));
% SetUp a whole lot of indices and needed Matrices
% Reactions which are allowed to vary freely: missing DeltaG_r std
    idx.deltaG_allowedToVary=idx.noDeltaE_Rxns_noLumped;
% Reactions with known DeltaG_r std (without Exchange rxns, as their
% deltaG_std is per defintion not defined)
    idx.known=setdiff(1:n.Rxns,idx.noDeltaE_Rxns);
    n.known=length(idx.known);
    [ExtracellRxnList] = findExcRxns(prepped_model);
    idx.internal=find(~ExtracellRxnList)';
    n.internal=length(idx.internal);
    idx.knWOexc=intersect(setdiff(idx.known,idx.deltaG_allowedToVary),idx.internal);
    idx.known = idx.knWOexc;
% Get matrices for TMFA constraints
    N=[S(:,idx.known) L]; %L is the stoichiometric matrix for lumped rxns
    idx.lumped=[length(idx.known)+1:size(N,2)]';
    n.lumped=size(a,1);
    n.knlu=size(N,2); %knlu = all rxns with known DeltaG0 (incl. lumped)
    idx.knlu=[idx.known idx.lumped'];
% different identitiy matrices
    AllMat=eye(n.Rxns+n.lumped,n.Rxns+n.lumped);
    %only take the rows of rxns which DeltaG is known (original+lumped)
    matG_a=AllMat(idx.knlu,:);
    %only take the rows of Lumped reactions
    matG_l=AllMat(idx.lumped,:); 
    % all reactions
    I_r=eye(n.Rxns); 
    % all lumped reactions
    I_l=eye(n.lumped);
    % metabolites that allowed to violate their default boundaries
    I_e=eye(n.Species);
    essSpecies=I_e(idx.essentAssocMets,:);
%% Boundaries
    %Fluxes
    ub=prepped_model.ub;
    lb=prepped_model.lb;
    %DeltaG_r
    Gmin=repelem(-300,n.Rxns+n.lumped);%-K
    Gmax=repelem(300,n.Rxns+n.lumped);%K
    %binaries 
    bin_y_lb=zeros(n.lumped,1);
    bin_z_lb=zeros(n.Rxns,1);
    bin_y_ub=ones(n.lumped,1);
    bin_z_ub=ones(n.Rxns,1);

    %Fix ATP/ADP ratio to 30 % can be enabled - if the concentration is not
    %fixed
    FixRatio=zeros(n.Species,1);
    %FixRatio(findMetIDs(prepped_model,'atp[c]'))=1;
    %FixRatio(findMetIDs(prepped_model,'adp[c]'))=-30;
    
%% Equality Constraints
    % Aeq     v                            deltaGr                       ln(conc)                   binary_y(lumped)           binary_z(all)        alpha(essentialMets)      beta(essentialMets)                                                                         
    Aeq=[ S                    zeros(n.Species,n.Rxns+n.lumped)  zeros(n.Species,n.Species)  zeros(n.Species,n.lumped)  zeros(n.Species,n.Rxns)     zeros(n.Species,n.essen)  zeros(n.Species,n.essen) ;
         zeros(n.knlu,n.Rxns)  matG_a                            -R*T*N'                     zeros(n.knlu,n.lumped)     zeros(n.knlu,n.Rxns)        zeros(n.knlu,n.essen)     zeros(n.knlu,n.essen)    ;
         zeros(1,n.Rxns)       zeros(1,n.Rxns+n.lumped)          FixRatio'                   zeros(1,n.lumped)          zeros(1,n.Rxns)             zeros(1,n.essen)          zeros(1,n.essen)        ]; 
    %beq
    beq=[zeros(n.Species,1); prepped_model.DeltaG_r_total(idx.known);lumpedModel.DeltaG_r_total;0];
%% Inequality Constraints
    % adding optimality constraint as well as adjusting concentration
    % constrait with the introduction of soft boundaries
    %A     v                            deltaGr                       ln(conc)                   binary_y(lumped)           binary_z(all)       alpha(essentialMets)    beta(essentialMets) 
    A=[-I_r                      zeros(n.Rxns,n.Rxns+n.lumped)   zeros(n.Rxns,n.Species)    zeros(n.Rxns,n.lumped)    zeros(n.Rxns,n.Rxns)      zeros(n.Rxns,n.essen)   zeros(n.Rxns,n.essen);
        I_r                      zeros(n.Rxns,n.Rxns+n.lumped)   zeros(n.Rxns,n.Species)    zeros(n.Rxns,n.lumped)    diag(-ub)                 zeros(n.Rxns,n.essen)   zeros(n.Rxns,n.essen);
        I_r                      zeros(n.Rxns,n.Rxns+n.lumped)   zeros(n.Rxns,n.Species)    zeros(n.Rxns,n.lumped)    zeros(n.Rxns,n.Rxns)      zeros(n.Rxns,n.essen)   zeros(n.Rxns,n.essen); 
        zeros(n.knlu,n.Rxns)     matG_a                          zeros(n.knlu,n.Species)    zeros(n.knlu,n.lumped)    I_r(idx.knlu,:).*K        zeros(n.knlu,n.essen)   zeros(n.knlu,n.essen); 
        zeros(n.lumped,n.Rxns)   matG_l                          zeros(n.lumped,n.Species)  I_l.*K                    zeros(n.lumped,n.Rxns)    zeros(n.lumped,n.essen) zeros(n.lumped,n.essen);
        zeros(n.lumped,n.Rxns)   zeros(n.lumped,n.Rxns+n.lumped) zeros(n.lumped,n.Species)  -I_l                      a                         zeros(n.lumped,n.essen) zeros(n.lumped,n.essen);
        zeros(n.essen,n.Rxns)    zeros(n.essen,n.Rxns+n.lumped)  essSpecies                 zeros(n.essen,n.lumped)   zeros(n.essen,n.Rxns)     -eye(n.essen)           zeros(n.essen,n.essen);   %conc - alpha <= concUB
        zeros(n.essen,n.Rxns)    zeros(n.essen,n.Rxns+n.lumped)  -essSpecies                zeros(n.essen,n.lumped)   zeros(n.essen,n.Rxns)     zeros(n.essen,n.essen)  -eye(n.essen)         ;   %-conc - beta <= -concLB
        -obj'                    zeros(1,n.Rxns+n.lumped)        zeros(1,n.Species)         zeros(1,n.lumped)         zeros(1,n.Rxns)           zeros(1,n.essen)        zeros(1,n.essen)       ]; %optimality constraint
        
    %b
    epsilon=0;
    b=[zeros(2*n.Rxns,1);ub; ones(n.knlu+n.lumped,1).*K-epsilon ; sum(a,2)-1;prepped_model.conc_ub(idx.essentAssocMets);-prepped_model.conc_lb(idx.essentAssocMets);-optVal*percentage]; 

%% Put everything together
    essentAssocMet.lb=prepped_model.conc_lb;
    essentAssocMet.lb(idx.essentAssocMets)=-Inf;
    essentAssocMet.ub=prepped_model.conc_ub;
    essentAssocMet.ub(idx.essentAssocMets)=Inf;
    
    lb_all=[lb;Gmin'; essentAssocMet.lb ; bin_y_lb; bin_z_lb; zeros(n.essen*2,1)];                     
    ub_all=[ub;Gmax'; essentAssocMet.ub; bin_y_ub; bin_z_ub; ones(n.essen,1)*250;ones(n.essen,1)*250]; 
    
    a_all=vertcat(Aeq,A);
    b_all=vertcat(beq,b);
                                                 
    vartype=[repmat('C', n.Rxns, 1); repmat('C', n.Rxns+n.lumped, 1);repmat('C', n.Species, 1); repmat('B', n.lumped, 1); repmat('B', n.Rxns, 1); repmat('C', n.essen*2, 1) ];
    ctype = [repmat('S', 1, size(Aeq,1)) repmat('U', 1, size(A,1))];
%% new objective
    c= [zeros(n.Rxns+n.Rxns+n.lumped+n.Species+n.lumped+n.Rxns,1);ones(n.essen+n.essen,1)];

%% RUN TMFA
% intlinprog
    %% intlinprog
    intcon=(n.Rxns+n.Rxns+n.lumped+n.Species+1):(n.Rxns+n.Rxns+n.lumped+n.Species+n.lumped+n.Rxns);
    options = optimoptions('intlinprog','Display','iter','MaxTime',480);
    
    [X,FVAL,STAT] =intlinprog(c,intcon,A,b,Aeq,beq,lb_all,ub_all,options);
    disp('Intlinprog Used')
    FVAL
    
    if(~isempty(FVAL))
%% isolate the parameters alpha and beta
        alpha_para= X((n.Rxns+n.Rxns+n.lumped+n.Species+n.lumped+n.Rxns+1):end-n.essen);
        beta_para = X((n.Rxns+n.Rxns+n.lumped+n.Species+n.lumped+n.Rxns+n.essen+1):end);

        prepped_model.conc_lb(idx.essentAssocMets)=prepped_model.conc_lb(idx.essentAssocMets)-beta_para;
        metsConcRanges.lb=prepped_model.conc_lb;

        prepped_model.conc_ub(idx.essentAssocMets)=prepped_model.conc_ub(idx.essentAssocMets)+alpha_para;
        metsConcRanges.ub=prepped_model.conc_ub;

        changedMets=prepped_model.mets(idx.essentAssocMets(unique([find(alpha_para~=0) ; find(beta_para~=0)])));

        for changeMet=1:length(changedMets)
            fprintf('For Met %s the Boundaries are changed to %d : %d\n',cell2mat(changedMets(changeMet)),(metsConcRanges.lb(findMetIDs(prepped_model,changedMets(changeMet)))),(metsConcRanges.ub(findMetIDs(prepped_model,changedMets(changeMet)))))    
        end
%% inital solution to give to TMFA
        x0=X(1:(n.Rxns+n.Rxns+n.lumped+n.Species+n.lumped+n.Rxns));
    else
        metsConcRanges=[];
        x0=[];
        fprintf('Empty Concentration Ranges and empty initial solution provided\n')
    end
end