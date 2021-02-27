%% TMFA
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
%     optional: provide x0 as 4th argument
% OUTPUT:
%     X
%        output from intlinprog TMFA program
%     FVAL
%        output from intlinprog TMFA program
%     STAT
%        output from intlinprog TMFA program
%     usedModel
%        constraint Matrix
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, FVAL, STAT,usedModel]=TMFA_LS(varargin)
%% Parameter & Input Settings 
prepped_model=varargin{1};
bigM=varargin{2};
idxtoBelumped=varargin{3};

if length(varargin)>3
    x0=varargin{4};
else
    x0=[];
end

original=prepped_model;
[n.Species,n.Rxns]=size(prepped_model.S);
R=constant.R_kcal;          
T=constant.T0;                  %[K] Temperature given
K=bigM;                        %upper boundary on DeltaG
deltaPsi=constant.deltaPsi;    %[mV]
deltapH=constant.deltapH;

%% Lumping
% search for metabolites with unknown Energy as they cause reactions to
% have undetermined deltaG values; can be skipped if idxtoBelumped is empty
  if ~isempty(idxtoBelumped)
        if iscell(prepped_model.DeltaG_m_std)
            prepped_model.DeltaG_m_std=cell2mat(prepped_model.DeltaG_m_std);
        end
        
        [alpha,lumpedModel,~]=lumpReactions(prepped_model,idxtoBelumped,1);
        
        % Update lumped Model with metNames + Energy
        lumpedModel.metNames=prepped_model.metNames;
        lumpedModel.DeltaG_m_std=prepped_model.DeltaG_m_std;
        if(isfield(prepped_model,'pH_offset'))
            lumpedModel.pH_offset=prepped_model.pH_offset;
        end
        if(isfield(prepped_model,'metCharges'))
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
       
    % Assinging deltaG Lumping    
    if ~isempty(L)
        disp('Assinging deltaG to Lumped Rxns')
        rxnIDXtoAssign=1:size(lumpedModel.S,2);
        lumpedModel.unKnownG_m=isnan(lumpedModel.DeltaG_m_std);
        %init result vector all NaN (for lumping reactions no NaN should be left after running the next for loop!)
        lumpedModel.DeltaG_r_std=nan(n.lumped,1);
        lumpedModel=calcDeltaG_r(lumpedModel);   
    else
         lumpedModel.DeltaG_r_total=[];
    end
    %adjust with tolerance
    tol=1e-10;
    lumpedModel.DeltaG_r_total(abs(lumpedModel.DeltaG_r_total)<tol)=0;
    
%% Find out reactions that are not constraint and not included in any lumped reactions 
% this reactions are allowed to vary freely with respect to their deltaG
    idx.noDeltaE_Rxns=find(isnan(prepped_model.DeltaG_r_total));
    idx.noDeltaE_Rxns=[idx.noDeltaE_Rxns; find(prepped_model.c==1)]; %biomass deltaG is not calculated
        
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
    obj=prepped_model.c;
    disp('Objective Function to be maximized is:');
    idx.Bio=find(prepped_model.c==1);
    disp(prepped_model.rxns(idx.Bio));
    
    %% SetUp a whole lot of indices and needed Matrices
    %Reactions vals which are allowed to vary freely + (Biomass)
    idx.deltaG_allowedToVary=idx.noDeltaE_Rxns_noLumped;
    [n.Species,n.Rxns]=size(S);
    idx.known=setdiff(1:n.Rxns,idx.noDeltaE_Rxns);
    n.known=length(idx.known);
    
    [ExtracellRxnList] = findExcRxns(prepped_model);%Delta G only for internal 
   
    
    idx.internal=find(~ExtracellRxnList)';
    n.internal=length(idx.internal);
    idx.knWOexc=intersect(setdiff(idx.known,idx.deltaG_allowedToVary),idx.internal); %Delta G known without exchange reactions

    idx.known = idx.knWOexc;
    
    %get reactions that get constraint
    N=[S(:,idx.known) L]; %L is the stoichiometric matrix for lumped rxns
    idx.lumped=[length(idx.known)+1:size(N,2)]';
    n.lumped=size(a,1);
    n.knlu=size(N,2); %knlu = all rxns with known DeltaG0 (incl. lumped)
    idx.knlu=[idx.known idx.lumped'];
    
    % Helping matrix
    AllMat=eye(n.Rxns+n.lumped,n.Rxns+n.lumped);
    matG_a=AllMat(idx.knlu,:);%only take the rows of rxns which DeltaG is known (original+lumped)
    matG_l=AllMat(idx.lumped,:); %only take the rows of Lumped reactions
    
    %IdentityMatrices
    I_r=eye(n.Rxns);
    I_l=eye(n.lumped);
    %% Boundaries
    %Fluxes
    ub=prepped_model.ub;
    lb=prepped_model.lb;
    
    %DeltaG_r
    Gmin=repelem(-300,n.Rxns+n.lumped);%must be greater than -K
    Gmax=repelem(300,n.Rxns+n.lumped);%must be smaller than K
    %binaries 
    bin_y_lb=zeros(n.lumped,1);
    bin_z_lb=zeros(n.Rxns,1);
    bin_y_ub=ones(n.lumped,1);
    bin_z_ub=ones(n.Rxns,1);
%% Equality Constraints
    % Aeq     v                            deltaGr                       ln(conc)                   binary_y(lumped)           binary_z(all)                                                                              
    Aeq=[ S                    zeros(n.Species,n.Rxns+n.lumped)  zeros(n.Species,n.Species)  zeros(n.Species,n.lumped)  zeros(n.Species,n.Rxns);
         zeros(n.knlu,n.Rxns)  matG_a                            -R*T*N'                     zeros(n.knlu,n.lumped)     zeros(n.knlu,n.Rxns)   ];
    %beq
    beq=[zeros(n.Species,1); prepped_model.DeltaG_r_total(idx.known);lumpedModel.DeltaG_r_total];

%% Inequality Constraints
    %A     v                            deltaGr                       ln(conc)                   binary_y(lumped)           binary_z(all)                                            %%
    A=[-I_r                      zeros(n.Rxns,n.Rxns+n.lumped)   zeros(n.Rxns,n.Species)    zeros(n.Rxns,n.lumped)    zeros(n.Rxns,n.Rxns);
        I_r                      zeros(n.Rxns,n.Rxns+n.lumped)   zeros(n.Rxns,n.Species)    zeros(n.Rxns,n.lumped)    diag(-ub);
        I_r                      zeros(n.Rxns,n.Rxns+n.lumped)   zeros(n.Rxns,n.Species)    zeros(n.Rxns,n.lumped)    zeros(n.Rxns,n.Rxns); 
        zeros(n.knlu,n.Rxns)     matG_a                          zeros(n.knlu,n.Species)    zeros(n.knlu,n.lumped)    I_r(idx.knlu,:).*K ; 
        zeros(n.lumped,n.Rxns)   matG_l                          zeros(n.lumped,n.Species)  I_l.*K                    zeros(n.lumped,n.Rxns);
        zeros(n.lumped,n.Rxns)   zeros(n.lumped,n.Rxns+n.lumped) zeros(n.lumped,n.Species)  -I_l                      a                    ];
    %b
    epsilon=0;
    b=[zeros(2*n.Rxns,1);ub; ones(n.knlu+n.lumped,1).*K-epsilon ; sum(a,2)-1];

%% Put everything together
    lb_all=[lb;Gmin'; prepped_model.conc_lb; bin_y_lb; bin_z_lb];
    ub_all=[ub;Gmax'; prepped_model.conc_ub; bin_y_ub; bin_z_ub];
    
    a_all=vertcat(Aeq,A);
    b_all=vertcat(beq,b);
    %                                             
    vartype=[repmat('C', n.Rxns, 1); repmat('C', n.Rxns+n.lumped, 1);repmat('C', n.Species, 1); repmat('B', n.lumped, 1); repmat('B', n.Rxns, 1)];
    ctype = [repmat('S', 1, size(Aeq,1)) repmat('U', 1, size(A,1))];
    %
    c=[obj;zeros(n.Rxns+n.lumped+n.Species+n.lumped+n.Rxns,1)];
    c=c.*-1; %we maximixe
    
%% RUN TMFA
    %either intlinprog or glpk
    %% intlinprog
    intcon=(n.Rxns+n.Rxns+n.lumped+n.Species+1):size(c,1);
    options = optimoptions('intlinprog','Display','iter','MaxTime',86400);%86400 Time limit of 1 day per calc
    %options = optimoptions('intlinprog','Display','iter','LPPreprocess','none');
    [X,FVAL,STAT] =intlinprog(c,intcon,A,b,Aeq,beq,lb_all,ub_all,x0,options);
    disp('Intlinprog Used')
    FVAL
    
    %% glpk
    %param.itlim=100;
    %param.tmlim=120;
    %s=-1; % objective already .*-(1)%   
    %[X, FVAL, STAT]=glpk(c,a_all,b_all,lb_all,ub_all,ctype,vartype,s);
    %disp('GLPK Used')
    %FVAL
%% Additional OUTPUT
    usedModel.Aeq=Aeq;
    usedModel.beq=beq;
    usedModel.A=A;
    usedModel.b=b;
    usedModel.N=N;
    usedModel.constraintsLHS=a_all;
    usedModel.constraintsRHS=b_all;
    usedModel.vartype=vartype;
    usedModel.ctype=ctype;
    usedModel.lb=lb_all;
    usedModel.ub=ub_all;
    usedModel.obj=c;
    usedModel.n=n;
    usedModel.idx=idx;
    usedModel.intcon=intcon;
    usedModel.K=K;
end