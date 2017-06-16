function Params = setup_params

% Setup the life history and simulation parameters for the Caribbean
% lobster model

% (I have just put in dummy numbers as placeholders)

% Age: (this is an age-structured modeled)
Params.Maxage = 16; % this can be an arbitrarily large number, Kanciruk 1980
Params.A = 1:Params.Maxage; % age vector 

% Growth: assume a von Bertalanffy function
Params.k = 0.24; % units = year^-1 - ok, Leon et al 2005
Params.Linf = 183.55; % asymptotic max length (units = cm or mm) - in mm, Leon et al 2005
Params.t0 = 0.45; % effective age at which L = 0, Leon et al 2005
%Params.L = Params.Linf .* (1 - exp(-k.*(Params.A - Params.t0))); 
Params.L = Params.Linf .* (1 - exp(-Params.k.*(Params.A - Params.t0))); %I modified k by the value of Params.k, check if ok
%Params.LW = 0.1 .* Params.L .^ 3; % convert length into biomass
Params.LW = [0.381,0.78,1.097,1.349,1.575,1.745,1.904,2.050,2.19,2.296,2.395,2.483,2.583,2.655,2.732*ones(1,Params.Maxage-14)];
% Based on SEDAR estimates from L-W sampling of tagged lobsters

% Fishing fleet 
Params.fleet_info = 1; % values range from 1 to infinity. Larger numbers = fishermen have less info, fishing effort is spread more evenly along coast

% Mortality (natural):
Params.M  = 0.34; % units = year^-1 - ok (Delury model base run, SEDAR 2010)
%M=Params.M; %I added that
% Mortality (harvest):
Params.Lf = 3; % age at which they enter the fishery - ok SEDAR 2010

% Maturity/reproduction
%Params.Amat = 2.5; %minimum age at maturity - ok 2-3 years (from book Lobster edited by Bruce F. Phillips)
%Params.Fec = 1 .* Params.LW .* (Params.A >= Params.Amat); % convert biomass into # eggs produced, for mature individuals only
% All of the following from SEDAR:
Mat = [0 0.5 0.75 ones(1,Params.Maxage-3)];
Broods = [ones(1,6), 2*ones(1,Params.Maxage-6)];
Eggs = [70983,270292,406946,507239,592776,654809,711260,761798,809349,844449,877102,905947,938062,960976,985381*ones(1,Params.Maxage-14)]; 
Params.Fec = Mat.*Broods.*Eggs;

% Density-dependent post-settlement survival:
% (assume a Beverton-Holt relationship)
Params.BH_alpha = 0.1; % survival at low density. R'(0)=0.1 creates a population where F20% is approximatly MSY (SEDAR, page 111 + equations White (2010))
Params.BH_beta = 2/500; % Maximum density per unit habitat (per m2)
% Based on max density observed in FWC 500 m^2 belt transects
% the alpha parameter could either be estimated from reef-scale
% experimental data, if such exist. Or we can scale-down a large-scale
% stock-recruit steepness parameter, if that exists. This will be a point
% of discussion for us.


% Calculate LEP: lifetime eggs produced per recruit:
LEP = sum(exp(-Params.M.*(Params.A-1)).*Params.Fec);

% Calculate FLEP: how LEP is reduced (proportionately) by harvest:
Fs = 0:0.01:2; % assume a harvest rate of 2 year^-1 is super high
FLEP = zeros(size(Fs));
for f = 1:length(Fs) %was length(F) before, check if the edit is correct 
    LEP1 = (exp(-Params.M.*(Params.A(Params.A<Params.Lf)-1)).*Params.Fec(Params.A<Params.Lf));
    LEP2 = exp(-Params.M*(Params.Lf-2)).*(exp(-(Params.M+Fs(f)).*(1:Params.Maxage-(Params.Lf-1))).*Params.Fec(Params.A>=Params.Lf));
    FLEP(f) = sum(LEP1) + sum(LEP2);
end

FLEP = FLEP./FLEP(1); % scale relative to  unfished

Params.LEP = LEP;
Params.Fs = Fs;
Params.FLEP = FLEP;





