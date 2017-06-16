function Lysel_test

% Try to figure out why deterministic & stochastic results are so far off
% for Lysel

% This file sort of mimics runme_lobster. Quick & dirty

% This code is not really used in the model per se, but I did write some
% snippets down below that are handy for plotting, could re-use if desired)

Type =1 ;

% Load stuff in:
Params = setup_params; 
Connmat = load_connmat(Type);
% Load habitat:
load('vector_reef_polygon.txt');
Hab=vector_reef_polygon;
%Nation=textread('vector_nation.txt','%s');
Nation = importdata('nation.csv');
Nation = Nation.data;
% Strip out 'fishing grounds' 
% This is only a guess, just taking out the appropriate number of cells
% from Honduras
Hon = find(Nation==4); % Honduras;
Hon2 = Hon(1:133);
Hon = Nation==4;
Hon(Hon2) = false;
Hab = Hab(~Hon);
Nation = Nation(~Hon);

% Plotting order:
PlotOrder = [37, 8, 7, 6, 5, 4, 2, 1, 9,...
             10, 11, 12, 3, 14, 15, 16, 17,...
             28, 18, 25, 27, 23, 27, 21,...
             20, 34, 19, 32, 33, 29, 30,...
             31, 24, 35, 36,38];

if size(Connmat,3) == 1
Params.eig = max(eig(Connmat));
else
Params.eig = max(eig(mean(Connmat,3))); 
end

Nats = length(unique(Nation));
CN = nan(Nats+1,Nats+1);
for n = 1:Nats
    for nn = 1:Nats
    OKn1 = Nation==n;
    OKn2 = Nation==nn;
    C = Connmat(OKn1,OKn2,:);
    CN(n,nn) = mean(C(:));
    end
end

% Re-create Figure 2A: (dispersal only)
CN(CN==0) = NaN;
pcolor(CN(PlotOrder,PlotOrder)); axis square

[SP,CM]= calculate_persistence(Params,Connmat,Nation,0.32*ones(length(Nation),1),Hab);

CN2 = nan(Nats+1,Nats+1);
for n = 1:Nats
    for nn = 1:Nats
    OKn1 = Nation==n;
    OKn2 = Nation==nn;
    C = mean(CM(OKn1,OKn2,:),3);
    CN2(n,nn) = mean(mean(C));
    end
end
figure
CN2(CN2==0) = NaN;
pcolor(CN2(PlotOrder,PlotOrder)); axis square; colormap jet

Connmat = load_connmat(2);
F = 0.34;
Fvec = F*ones(2930,1);
[N,Neq,CM,CMn,MP,SP,SR,LR,NP] = run_model(Params, Connmat, min(1,Hab), Nation, Fvec);

CN2 = nan(Nats+1,Nats+1);
for n = 1:Nats
    for nn = 1:Nats
    OKn1 = Nation==n;
    OKn2 = Nation==nn;
    C = mean(CM(OKn1,OKn2,:),3);
    CN2(n,nn) = mean(mean(C));
    end
end
figure
CN2(CN2==0) = NaN;
pcolor(CN2(PlotOrder,PlotOrder)); axis square; colormap jet

keyboard
    