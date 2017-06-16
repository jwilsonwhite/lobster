function [SP, CM, NP] = calculate_persistence(Params,Connmat,Nation,Fvec,Hab,Subnet)

% Calculate static self-persistence on a nation-by-nation basis

% Output: Vector SP with self-persistence for each nation (values >= 1
%               indicate self-persistence)
%         Matrix CM with connectivity to and from each nation

% Could also calculate network persistence for particular subnetworks, but
% would need to specify them somehow

% 1) Assemble projection matrix:

% FLEP for each nation:
FLEP = nan(length(Fvec),1);
for i = 1:length(Fvec)
    F_OK = abs(Fvec(i) - Params.Fs); % find the correct value of F
    F_OK = find(F_OK == min(F_OK),1);
    FLEP(i) = Params.FLEP(F_OK);
end

% Assemble the connectivity matrix:

% Use the following if BH_alpha is already scaled to LEP:
% Add correction to alpha value to account for slope of BH (White 2010 Fish
% Res)

Connmat = bsxfun(@times,Connmat,Hab'); % multiply columns by Hab


Alpha = 1./Params.BH_alpha./max(eig(Connmat)); % rescale slope (including Habitat effects)

CM = bsxfun(@times,Connmat,FLEP(:)')*Alpha;

% OR if BH_alpha has not been rescaled, use this:
%CM = Connmat .* repmat(FLEP(:),[1,length(FLEP)]) .*Params.LEP .* Params.BH_alpha;

% 2) Self-persistence 
% Loop over each nation, obtain the submatrix & calculate persistence:
Ns = length(unique(Nation)); % number of nations
SP = nan(Ns,1); % pre-allocate
for n = 1:Ns
    OKn = Nation==n; 
    CMsub = CM(OKn,OKn);
    SP(n) = max(eig(CMsub));
end

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


