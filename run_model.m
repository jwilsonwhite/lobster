function [N,Neq,CM,CMn,MP,SP,SR,LR,NP,BC] = run_model(Params,Connmat,Hab,Nation,Fvec,Subnet)

% Lobster population model

% Initialize the population vector
T = 500; % length of simulation. The model will be run for Tdd years to reach equilibrium, then 50 years w/o density dependence to test persistence (if necessary)
Tdd = 450; % time at which density-dependence turns off
P = length(Hab);

N = zeros(Params.Maxage,P,T);
N0 = zeros(Params.Maxage,P,Tdd); % to find the initial unfished conditions
Ft = zeros(P,T); % fishing at each time in each patch
E = zeros(P,T); % egg production at each time in each patch
uniNats = unique(Nation); % list of unique individual nations
Nats = length(unique(Nation)); % how many nations
whichC = nan(T,1); % track which connectivity matrix is used in each year

N0(1,:,1) = Params.BH_beta; % start off with a large settlement pulse

Alpha = 1/(Params.BH_alpha*Params.LEP)./Params.eig;


% Begin by finding initial (unfished) conditions, starting from arbitrary
% point:

% survival matrix:
A = eye(Params.Maxage-1)*exp(-Params.M);
A = [zeros(1,Params.Maxage-1); A];
A = [A,zeros(Params.Maxage,1)];


for t = 2:Tdd
    
    % Advance ages
    N0(:,:,t) = A*N0(:,:,t-1);
    
    % reproduction:
    E(:,t) = (Params.Fec*N0(:,:,t-1))'; % Px1 vector of eggs
    
    % dispersal:
    if size(Connmat,3)>1
    C = mean(Connmat,3); % for this initial period, use deterministic mean % (:,:,randperm(3,1)); % choose which matrix to use
    else
    C = Connmat;
    end
    L = C*E(:,t); % Px1 vector of settlers in each patch
    
    % Settlement mortality: 
    S = Alpha*L./(1 + Alpha*L./(Params.BH_beta.*Hab));

    % Add in new settlers
    N0(1,:,t) = S;
end % end loop over t

% Now we have unfished equilibrium:
N(:,:,1) = N0(:,:,end);

% Now simulate dynamics including fishing & dynamic fleet:
Ft(:,1) = Fvec;
BC = ones(size(Fvec));
CPUE = ones(size(Fvec));

for t = 2:T

    if t == Tdd
        N(:,:,t-1) =  N(:,:,t-1)*1e-3;
    end
    
    % Loop over each patch, create transition matrix with dynamic fleet:
    for n = 1:Nats % loop over each nation
        Nok = Nation == uniNats(n); % which patches
        CPUEn = CPUE(Nok,t-1).^(1/Params.fleet_info);
        Ftot = sum(Ft(Nok,t-1));
        Fnew = Ftot*CPUEn/sum(CPUEn); % redistribute based on prior year's CPUE
        Ft(Nok,t) = Fnew;
    end % end loop over nations
    
    % Loop over patches
    for i = 1:P
        
        A(2:end,1:end-1) = diag(exp(-(Params.M + Ft(i,t)*(Params.A(1:end-1)>=(Params.Lf-1)))));
        
        N(:,i,t) = A*N(:,i,t-1);
        
        % Biomass caught:
        BC_tmp = (N((Params.Lf-1):end-1,i,t-1) - N((Params.Lf):end,i,t)).*Ft(i,t)/(Ft(i,t)+Params.M);
        BC(i,t) = max(0,Params.LW((Params.Lf-1):end-1)*BC_tmp);
        
        
    end % end loop over patches
    
    % reproduction:
    E(:,t) = (Params.Fec*N(:,:,t))'; % Px1 vector of eggs
    
    % dispersal:
    if size(Connmat,3)>1
    whichC(t) = randperm(3,1); % keep track of which matrix was used    
    C = Connmat(:,:,whichC(t)); % choose which matrix to use
    else
    C = Connmat;
    end
    L = C*E(:,t); % Px1 vector of settlers in each patch
    
    % Settlement mortality: 
    if t > Tdd && size(Connmat,3)>1
    S = Alpha*L; % no DD 
    else
    S = Alpha*L./(1 + Alpha*L./(Params.BH_beta.*Hab));
    end
    
    % Add in new settlers
    N(1,:,t) = S;
    
    % CPUE:
    CPUE(:,t) = max(eps,BC(:,t)./Ft(:,t));
end % end loop over t

% At equilibrium, check for overall persistence & self-persistence

Lag = 10;

%%%%% Do different calculations depending on whether it is a deterministic
%%%%% or stochastic simulation:

if size(Connmat,3)==1 % constant environment - can do standard calculations
C = Connmat;

% assemble connectivity matrix with FLEP estimated from equilibrium fishing
% rate. (deprecated, this is a bad way to do it)
%Flepvec = interp1(Params.Fs,Params.FLEP,min(Ft(:,end),max(Params.Fs)));
%CM = C .* repmat(Hab(:)',[P,1]) .* repmat(Flepvec(:),[1,P])' .* Alpha .*Params.LEP;
%CM(isnan(CM))=0;

% Alternative using Eggs directly:
EPR = E(:,Tdd-Lag)./squeeze(N(1,:,Tdd-Lag))'; % eggs per recruit
EPR(isnan(EPR)) =0;
EPR(isinf(EPR)) =0;
CM = C .* repmat(EPR(:)',[P,1]) .* Alpha;
CM(isnan(CM)) = 0;

% Persistence of entire metapopulation:
MP = max(eig(CM));

% Loop over each nation:
SP = nan(Nats,1);
SR = SP; LR = SP;
for n = 1:Nats
    OKn = Nation==n;
    CMsub = CM(OKn,OKn);

    SP(n) = max(eig(CMsub));  % self-persistence
    
    SR(n) = sum(diag(CMsub))./sum(CMsub(:)); % self-recruitment
    
    LR(n) = sum(diag(CMsub))./sum(sum(CM(:,OKn))); % local retention
    
end %end loop over nations

% 3) Network persistence (if the input argument Subnet was provided)
if exist('Subnet','var')
   NNs = length(Subnet); % number of subnetworks
   NP = nan(NNs,1);
   for n = 1:NNs
       Subn = Subnet{n}; % vector of subnetwork nations
       OKn = false(size(Nation));
       for s = 1:length(Subn)  % loop over each nation to figure out which nodes are in the subnetwork
           OKn = OKn | Nation == Subn(s);
       end
       CMsub = CM(OKn,OKn); % now it is just the same calculation as before
       NP(n) = max(eig(CMsub));
   end
else
    NP = NaN;
end



else % if we are using multiple connectivity matrices, then the eigenvector-based calculations are not possible.

    SP = NaN; NP = NaN; % these are no longer calculated, so put in dummy placeholders.
    
    % The whole connectivity matrix is too large to calculate the
    % covariances necessary to use Tuljapurkar's approximation for the
    % stochastic growth rate. So we have to use the simulation
    % approximation. This will be done with the final 50 years of the
    % simulation that was run without density-dependence
    
    if P == 1 % (some nations just have 1 patch)
    Nsum = squeeze(sum(N(:,:,Tdd:T),1));
    else
    Nsum = sum(squeeze(sum(N(:,:,Tdd:T),1)),1);
    end
 
    r = diff(log10(Nsum)); % this is equivalent to log lambda (log of the annual growth rate)
    
    MP = exp(mean(r)); % this is exp(loglambdaS), the log stochastic long term growth rate. 
                       % Exponentiated to get in the same units as the static eigenvector calculation [using Caswell's nomenclature]

if isnan(MP); MP=0; end % this happens if the population has declined to within the numerical tolerances of Matlab.

    % Now make the calculations of self-recruitment and local retention.
    % To do this, assemble a series of connectivity matrices with eggs-per-recruit. Do
    % this only for the part of the simulation before DD was removed
    Eggs = mean(E(:,Tdd-Lag-(1:Lag)-2),2);
    if size(Eggs)>1 % if there is more than one polygon in this nation
    Rec = mean(squeeze(N(1,:,Tdd-Lag-(1:Lag)-2)),2); % eggs per recruit
    else
    Rec = mean(N(1,1,Tdd-Lag-(1:Lag)-2));
    end
    EPR = Eggs./Rec;
    EPR(isnan(EPR)) =0;
    EPR(isinf(EPR)) =0;
    
    CM = nan(size(Connmat,1),size(Connmat,2),Lag);
    for l = 1:Lag
        Index = Tdd - (Lag-l+2);
        
        % Old stuff, now deprecated:
       % Flepvec = interp1(Params.Fs,Params.FLEP,min(Ft(:,Index),max(Params.Fs)));
       % CM(:,:,l) = Connmat(:,:,whichC(Index)) .* repmat(Hab(:)',[P,1])  .* repmat(Flepvec(:),[1,P])' .* Alpha .*Params.LEP; 
      %  EPR = E(:,Index)./squeeze(N(1,:,Index))'; % eggs per recruit
      %  EPR(isnan(EPR)) =0;
      %  EPR(isinf(EPR)) =0;
      
        CM(:,:,l) = Connmat(:,:,whichC(Index)) .* repmat(EPR(:)',[P,1]) .* Alpha;
      
    end % end loop over l
    CM(isnan(CM)) = 0;
    
    Cmean = mean(CM,3); % average matrix
       
    % Loop over each nation:
    SR = nan(Nats,1); LR = SR;
    for n = 1:Nats
    OKn = Nation==n;
    CMsub = Cmean(OKn,OKn);
    
    SR(n) = sum(diag(CMsub))./sum(CMsub(:)); % self-recruitment
    
    LR(n) = sum(diag(CMsub))./sum(sum(CM(:,OKn))); % local retention
 %   
    end %end loop over nations

end % end if a constant connectivity matrix


Neq = mean(N(:,:,Tdd:Lag),3)./N0(:,:,end); % equilibrium size, relative to unfished

% Sum up connectivity matrix by nations for further analysis:
CMn = nan(Nats);
for n = 1:Nats
    for m = 1:Nats
CMn(n,m) = sum(sum(CM(Nation==n,Nation==m)));
    end
end

