function runme_lobster(F,type)

% Wrapper code to run all the components of the Caribbean lobster model

% Input arguments:
% F: fishing mortality rate (spatial average). Target values: 0, 0.34, 1
% type: 1 = average matrix (averaged over 3 years; deterministic run); 2 = use individual
% years (stochastic run)

savename = 'lobster_scenario1_testingJun2017.mat';
% Adjust this to keep track of the runs you make. This file will store the
% parameters, matrices, and model results for a particular run.

Subnetworks = {[1,2,4],... % mesoamerica
               [10,12],... % N Cuba, Bahamas
               [12,13],... %Cuba
               [13, 14, 3],...; % S Cuba, Jamaica, Cayman
               [9,17],... % USA
               [18,23,25,34,35],... % netherlands
               [18:27],...  % leeward islands
               [28:33],...    % windward islands
               [18:33]}; % leeward + windward

Params = setup_params; 
% this code will create all of the necessary life history parameters.
% It will return a structure file (Params) that can be passed to other
% functions

%type = 1;
Connmat = load_connmat(type); % Type = 1 (1 average matrix) or 2 (3 separate years)
% Load in the connectivity matrix (may need to add in possibility of
% different alternatives here, e.g., overall mean vs. multiple years.

% In general Connmat should be a p x p matrix (where p is the number of
% habitat patches in the model). Each element Connmat(i,j) gives the
% probability of dispersal from patch j (columns) to patch i (rows)

load_habitat;
% Hab should be a vector containing the proportion (or total area) of
% suitable reef habitat in each patch. It will be a p x 1 vector.
% (Need to decide whether it is the proportion out of a maximum or just an
% area, depending on whether the polygons themselves vary in size)

% It will be desirable to also create a p x 1 vector 'Nation' that indicates which
% nation/jurisdiction/island, etc. each patch belongs to. This will
% facilitate calculating self-persistence and whatnot. Each Nation will
% have a numerical code (1, 2, 3, ...)

% For plotting purposes, you will also want to have a p x 1 cell array that
% contains the corresponding vertices (in lat lon coordinates). That will
% make it easy to plot results.

% Find the values of F to use (just look at a one-patch population)

%for f = 1:length(Params.Fs)
%    [~,~,~,~,~,~,~,~,BC] = run_model(Params, 1, 1, 1, Params.Fs(f));
%    BCx(f) = BC(end);
%end
%Params.Fs(BCx==max(BCx)); % MSY

% Based on the analysis above, I recommend using F = 0 (unfished), F = 0.34
% (F_MSY), and F = 1.0 (highly overexploited)


Fvec = F*ones(4921,1); %Fvec = F*ones(4921,1)
% Fvec will be a p x 1 vector giving the fishing rate in each patch. This
% should be an annual instantaneous rate (units year^-1). If you have
% something like proportion of stock harvested, then we can discuss how to
% convert it into what we need. Presumably these values will vary on a
% nation-by-nation basis?

% Calculate the eigenvalue to save time later:
if size(Connmat,3) == 1
Params.eig = max(eig(Connmat));
else
Params.eig = max(eig(mean(Connmat,3))); 
end

[SP, CM, NP] = calculate_persistence(Params, Connmat, Nation, Fvec, Hab,Subnetworks);
% This code will take the variables loaded in above and makes long-term,
% deterministic calculations about self-persistence (which nations are
% self-persistent), network persistence, and source-sink linkages, etc.
% These will be static calculations that assume no variability in
% connectivity, fishing, etc. 
% Input 'Subnetworks' is a cell array containing vectors listing
% subnetworks (nations) for which we want to see if they are a persistent
% subnetwork (it is infeasible to calculate this for all possible
% subnetworks, so we have to be targeted) (this is an optional entry)
% The output SP will be a structure variable with the various results.
% The output CM will be a NxN matrix, where N is the number of nations,
% giving the connectivity (# offspring per generation) among each patch
% (this will be useful in looking at source-sink relationships)
% The output NP gives the network persistence for each set of nations in
% Subnetworks

if type == 1
% If we are using a single deterministic Connmat, we can make all of the
% calculations in one step using the eigenvalue approach.
[N,Neq,CM,CMn,MP,SP,SR,LR,NP] = run_model(Params, Connmat, Hab, Nation, Fvec, Subnetworks);

% This code will run a dynamic model. This can handle year-to-year
% variation in connectivity, as well as a dynamic fishing fleet (i.e., the
% fishing effort can move around within a jurisdiction to concentrate on 
% the more productive patches).
% The model will make calculations about self persistence, network
% persistence, source-sink relationships, and long-term patterns of biomass
% and yield.
% Outputs:
% N: age x patch x time vector of lobster density
% Neq: N, scaled to unfished equilibrium
% CM2: Connectivity matrix (dispersal matrix * egg production * slope of
% Bev-Holt)
% CMn: Connectivity matrix, aggregated at national scale
% MP: Persistence of entire metapopulation
% SP2: Self-persistence of each nation
% SR2: Self recruitment of each nation (proportion of settlers that were
% locally spawned)
% LR2: Local recruitment (proportion of larvae retained) for each nation
% NP2: Network persistence (if Subnetworks are provided)


elseif type ==2 % If multiple Connmats are used, the model is stochastic and an alternative set of calculations are necessary
    
% Overall simulation:    (this takes a while, so I recommend using %% to comment this out if you want to just look at self-persistent patterns) 
[N,Neq,CM,CMn,MP,~,SR,LR,~,BC] = run_model(Params, Connmat, Hab, Nation, Fvec, Subnetworks);    
    
% With a stochastic matrix, we cannot calculate persistence directly from
% the eigenvalues, because they vary from year to year. Instead, estimate long-term 
% growth rate (without density-dependence) by simulation
Nats = length(unique(Nation));

for n = 1:Nats
    OKn = Nation==n;
    for j = 1:10
    [~,~,~,~,SP2(n,j)] = run_model(Params, Connmat(OKn,OKn,:), Hab(OKn), Nation(OKn), Fvec(OKn));
    end
    
    C = Connmat(OKn,OKn,:);
    for i = 1:size(C,3)
        Eig(i) = max(eig(C(:,:,i)));
    end
    
        
    % Calculate CV of eigenvalues: (could also calculate some different
    % stats here if desired)
    % (look for correspondence between CV of eigs & difference in
    % deterministic vs. stochastic outcomes)
    CV_c(n) = std(Eig)./mean(Eig);


end
    

% Persistence calculation for individual subnetworks:
for n = 1:length(Subnetworks)
   Subn = Subnetworks{n};
   OKn = false(size(Nation));
   for s = 1:length(Subn)
       OKn = OKn | Nation == Subn(s);
   end
    [~,~,~,~,NP(n)] = run_model(Params, Connmat(OKn,OKn,:), Hab(OKn), Nation(OKn), Fvec(OKn));


end
end % end if Type 1 or 2


% 

save(savename) % Save all the results

keyboard % can delete this if you don't want to enter debug mode here.

% To plot the connectivity matrix:
CMn_tmp = CMn;
% Add NaNs to the edges so that all values are plotted
CMn_tmp(:,end+1) = NaN; CMn_tmp(end+1,:)=NaN;

% Make the plot:
%figure(1)
%clf
%pcolor(CMn_tmp); % checkerboard plot

% Some plotting options:
%set(gca,'xtick',0.5:1:40.5,'ytick',0.5:1:40.5); % add ticks
%set(gca,'tickdir','out','ticklength',[0.02 0.02]) % make ticks face out
%ylabel('Destination','fontsize',14); 
%xlabel('Origin','fontsize',14)
%ch = colorbar; % add colorbar
%set(ch,'tickdir','out','ticklength',[0.01 0.01])





