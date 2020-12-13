function [lumpVector,lumpedRxn,reactionsCoeff]=basicLump(varargin)
    Matrix=varargin{1};
    MetToZero=varargin{2};
    if(~isempty(varargin{3}))
        intCutMat=varargin{3};
    else
        intCutMat=[];
    end
    unknownIDX=MetToZero;
    Lumped=[];
    alpha={};
    S=Matrix;
    reactionsCoeff=[];
    metsToCurateManually=[];
    rxnsNames=[];
    i=MetToZero;  
    [m,n]=size(S);
    metsToLump=i;
    
    %use to restirct x vars of those rxns
    currentInvolvedRxns=find(sum(abs(S(metsToLump,:)),1)~=0);
    ub_x=zeros(n,1);
    ub_x(currentInvolvedRxns)=1000;
    ub=repmat(1000,m,1);
    lb_x=zeros(n,1);
    lb_x(currentInvolvedRxns)=1;
    
    Aeq=[S -eye(m) eye(m) zeros(m)];
    beq=zeros(m,1);
    A=[zeros(1,n) zeros(1,m) zeros(1,m) zeros(1,m)];
    b=[0];
    ub_all=[ub_x;ub;ub;ones(m,1)];
    lb=[lb_x;zeros(m,1);zeros(m,1);zeros(m,1)];
    obj=[zeros(1,n) ones(1,m) ones(1,m) zeros(1,m) ]';
    
    obj=obj.*1 ; % min
    Aeq_tmp=Aeq;
    beq_tmp=beq;
    % expand equality matrix Aeq for each metabolite of dependent
    % metabolite group
    if size(metsToLump,1)>1
        metsToLump=metsToLump';
    end
    for k=metsToLump
        cons_plus=zeros(1,m);
        cons_plus(k)=1;
        cons_minus=zeros(1,m);
        cons_minus(k)=-1;
        toAddAeq=[zeros(1,n) cons_plus  cons_minus zeros(1,m)];
        Aeq_tmp=[Aeq_tmp;toAddAeq];
        beq_tmp=[beq_tmp;0];
    end
    intcon=n+2*m+1:length(ub_all);
    
    if ~isempty(intCutMat)
        % Introduce Integer Cuts
        % n additional binary variables
        % Change of A, Aeq column dimension
        % Change of A and b row dimension (additional constraints)
        % Change objective dimension
        % Change lb and ub dimension
        % Change intcon
        Aeq_tmp=[Aeq_tmp zeros(size(Aeq_tmp,1),n)];
        obj=[obj; zeros(n,1)];
        lb=[lb; zeros(n,1)];
        ub_all=[ub_all; ones(n,1)];
        intcon=[intcon intcon(end)+1:intcon(end)+n];        
        
        addToA=[zeros(size(intCutMat,1),size(A,2)) intCutMat;
                zeros(size(intCutMat,1),size(A,2)) -intCutMat];
        A_tmp=A;
        A=[A_tmp zeros(size(A_tmp,1),n)]; 
        A=[A;addToA];
        b=[b;sum(intCutMat,2)-1;sum(intCutMat,2)+1];
    end
    
     options = optimoptions('intlinprog','MaxTime',120,'Display','off','IntegerPreprocess','none','IntegerPreprocess', 'none','LPPreprocess','none');
     [X,optVal,STAT]=intlinprog(obj,intcon,A,b,Aeq_tmp,beq_tmp,lb,ub_all,options);
        
     
     if STAT>0 
         %which rxns
         reactionsToLumped=logical(X(1:n)~=0); % get Support
         %get actual rxns coeffiecients
         reactionsCoeff=X(1:n);
         %Get new Lumped Reaction:
         tmp_newLumped=X(n+1:n+m)-X(n+m+1:n+m+m);
         if(nnz(tmp_newLumped(MetToZero,:))>0)
            disp('ISSUE!!!')
         end
         tmp_newLumped(abs(tmp_newLumped)<1e-6)=0;
         lumpedRxn=tmp_newLumped;
         lumpVector=[reactionsToLumped~=0];
     else
         %disp('No Solution Found')
         lumpedRxn=[];
         lumpVector={};
     end
end