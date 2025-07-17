% =====================================================================================================================================
% This script serves as the main controller for running stochastic simulations of Influenza A Virus (IAV)
% intracellular replication using the Gillespie algorithm with the tau-leaping approximation. It sets up the necessary parameters, initial conditions,
% and simulation settings, and then calls the script mode_IAV_GillespieTauLeaping. 
% =====================================================================================================================================

% =====================================================================================================================================
% This script is used to generate simulated data for the Figure 2B and Figure S1
% of the paper Segredo-Otero and Gresham 2025.
% =====================================================================================================================================

% The reaction network and all required parameters are defined externally in an Excel file, referred to as 'Param_mat'.
% This file must contain at least six sheets, each providing specific information needed for the simulation:

% 1. Stoichiometric matrix: A matrix that encodes how each reaction changes the copy number of molecular species.
%    Each row represents a species, and each column represents a reaction. Entries indicate the net gain or loss
%    of each species in that reaction.

% 2. Reactants matrix: Specifies which molecular species are reactants in each reaction. This is used to compute 
%    propensities (reaction probabilities per unit time) based on current species counts.

% 3. Rate constants matrix: A numerical matrix where rows correspond to individual kinetic parameters (e.g., k_on, k_off),
%    and columns correspond to reactions. Used to evaluate propensity functions for the stochastic simulation.

% 4. Parameter list: A list associating each kinetic parameter symbol with its corresponding value. This serves
%    both for readability and for parameter mapping when building the propensity functions.

% 5. QSSA-eliminated reactions: This sheet contains the parameters of reactions that have been removed from the 
%    Gillespie simulation under the Quasi-Steady-State Approximation (QSSA), i.e., fast reactions assumed to be in equilibrium.

% 6. Virion composition: Defines the number of each protein type (e.g., HA, NA, NP) found in a mature IAV virion. 
%    This is used to initialize the state of the system or to constrain assembly steps.

% Note: Individual degradation reactions for each species do not need to be explicitly included in the stoichiometric,
% reactants, or rate constants matrices. These degradation reactions are automatically generated within the script,
% but the corresponding degradation rate parameters must be defined in Sheet 4 (the parameter list).
% Upon execution, the script prompts the user to input two key parameters:

% - MOI (Multiplicity of Infection): Defined here as the number of genetically distinct IAV virions infecting the cell.
% - IniV (Initial Virions per Variant): Specifies the number of virions of each genetic variant that initiate infection.

% The total number of infecting virions will thus be MOI × IniV.

clear all
close all
clc

% ======================
% PARAMETERS AND INPUTS
% ======================

% Excel file containing the reaction network and all parameter definitions
Param_mat = 'Stock_matrix_MDCK.xlsx';

% Load the stoichiometric matrix (Sheet 1): rows = species, columns = reactions
MR = xlsread(Param_mat, 1);
sM = size(MR);  % Dimensions of the stoichiometric matrix

% Load the reactants matrix (Sheet 2): rows = species, columns = reactions
MS = xlsread(Param_mat, 2);

% Load the rate constants matrix (Sheet 3): rows = parameters, columns = reactions
MP = xlsread(Param_mat, 3);

% Load virion composition (Sheet 6): number of proteins per type in one IAV virion
NPs = xlsread(Param_mat, 6)';  % Transposed to be a row vector

% Load parameter list (Sheet 4): symbols and numerical values of model parameters
Params = xlsread(Param_mat, 4);

% Load QSSA-eliminated reaction parameters (Sheet 5)
Reg = xlsread(Param_mat, 5);

% ======================
% INITIAL CONDITIONS
% ======================

% Multiplicity of infection (MOI): number of genetically distinct virions infecting the cell
MOI = 10;

% Initial number of virions per genetic variant infecting the cell
IniV = 1;

% ======================
% SIMULATION SETTINGS
% ======================

% Frequency (in number of steps) at which the simulation outputs intermediate results
printFreq = 5e3;

% Maximum time step allowed between events (in minutes)
dtLim = 1;

% Maximum number of assembled virions allowed in simulation before termination
maxVir = 1e5;

% Threshold number of virions after which their count is stored for analysis
maxVir_store = 1e4;

% Maximum number of Gillespie steps allowed before forced termination
maxStep = 1e7;

% Maximum number of replication cycles stored
maxCycle = 30;

% Tolerance value for tau-leaping dt calculation
e = 0.03;

% Final simulation time (in minutes)
tf = 12*60;

% Total number of molecular species (from stoichiometric matrix)
NumS = sM(1,1);

% Initial molecular state vector (1 x NumS), initialized to zero
Si = zeros(1,NumS);

% Initial states for each of the MOI virions (MOI x NumS)
SiM = zeros(MOI,NumS);

% Extended matrix to store species states across multiple cycles (MOI x NumS * maxCycle)
SiMC = zeros(MOI,NumS*maxCycle);

% Initial simulation time
ti = 0;

% Number of IAV genomic segments (fixed at 8)
NumG = 8;

% Initial presence (1) of each segment in each virion (MOI x NumG)
V0 = ones(MOI, NumG);

% Probability of packing certain numbe of virions. Based on experimental
% data, 80% of IAV irions contain all segments, while 20% contain 7.
Frac_I = zeros(1,8); 
Frac_I(1) = 0.8;  % 80% of infection comes from variant 1
Frac_I(2) = 0.2;  % 20% from variant 2

% Minimum time resolution for simulation events (in minutes)
dtMin = 0.001;

%%
% Run the stochastic IAV replication model using a hybrid Gillespie/Tau-Leaping scheme
tic
[SiT, Tv, SiTM, SiMC, Virs_G, Virs_C, Virions_time, ...
 NumS, step, NumG, mRNA, Prots, nRNP, pRNP, PM_cyt, ...
 ti, br] = model_IAV_GillespieTauLeaping(MOI, IniV, tf, ...
 printFreq, maxVir, maxVir_store, maxStep, ...
 maxCycle, e, Si, SiM, SiMC, ti, V0, Frac_I, ...
 Params, Reg, MR, MS, MP, NPs, dtMin);
toc

%%        
% Reshape SiTM to 3D: [species x time x cell] for easier indexing
SiTM1 = zeros(NumS, step, MOI);
for m = 1:MOI
    IniS = (m-1)*NumS + 1;
    FinS = m*NumS;
    SiTM1(:,:,m) = SiTM(IniS:FinS,:);
end
SiTM = SiTM1;


% Sum of genomic RNA species across compartments (nRNP + pRNP + cytoplasmic)
TotalG = zeros(NumG, step, MOI);
Virs = zeros(1, step);
for m = 1:MOI
    TotalG(:,:,m) = SiTM(nRNP,:,m) + SiTM(pRNP,:,m) + SiTM(PM_cyt,:,m);
    Virs(1,:) = SiT(NumS,:);  % Total number of virions over time (last species in SiT)
end

%%
col = jet;
colM = floor(256/NumG);

figure
title('mRNA')
hold on
for g = 1:NumG
    plot(Tv,SiT(mRNA(g),:),'Color',col(g*colM,:),'Linewidth',2)
end
hold off
set(gca,'YScale','log')

figure
title('Genomes')
hold on
for g = 1:NumG
    plot(Tv,SiT(PM_cyt(g),:),'Color',col(g*colM,:),'Linewidth',2)
end
hold off
set(gca,'YScale','log')

    figure
    plot(Tv,SiT(Prots,:),'Linewidth',2)
    set(gca,'YScale','log')
title('Proteins')

colM = floor(256/MOI);
%%
figure
title('Genomes')
hold on
for m = 1:MOI
    plot(Tv,TotalG(:,:,m),'Color',col(m*colM,:),'Linewidth',2)
end
plot(Tv,sum(TotalG,3),'k','Linewidth',2)
hold off
set(gca,'YScale','log')

figure
title('Virions')
hold on
plot(Tv,Virs,'-k','Linewidth',2)
hold off
set(gca,'YScale','log')

%%

% Compute degradation probability based on competition with polymerase and NP binding

dRNA  = Reg(1);                   % vRNA degradation rate
PolB  = Reg(3) .* SiT(36,:);      % polymerase binding to vRNA
dRNAP = Reg(2);                   % degradation of vRNA-polymerase complex
NPB   = Reg(4) .* SiT(38,:);      % NP binding to vRNA

% Event probabilities
P_dRNA  = 1 - exp(-dRNA);
P_PolB  = 1 - exp(-PolB);
P_dRNAP = 1 - exp(-dRNAP);
P_NPB   = 1 - exp(-NPB);

% Normalized probabilities for competition
PPol = P_dRNA ./ (P_dRNA + P_PolB);     % Probability that vRNA is not protected by polymerase
PNP  = P_dRNAP ./ (P_dRNAP + P_NPB);    % Probability that vRNA-Pol is not protected by NP

% Combined probability of RNA being unprotected and degraded
P = PPol .* PNP;


%%
colM = jet(MOI); % Generar colores únicos

for s = 1:NumG
    figure % Crear una nueva figura en cada iteración
    title(['Segment ' num2str(s)])
    hold on
    
    % Mantener los mismos colores para cada línea en todas las figuras
    set(gca, 'ColorOrder', colM, 'NextPlot', 'replacechildren');
    
    plotG = squeeze(TotalG(s,:,:)); % Extraer la matriz correctamente
    plot(Tv, plotG, 'Linewidth', 2); % Graficar todas las líneas a la vez

    set(gca, 'YScale', 'log') % Escala logarítmica en Y
    grid on
    legend(arrayfun(@(x) sprintf('Virion %d', x), 1:MOI, 'UniformOutput', false),'NumColumns',2,'Location','SE');
    xlabel('Time (hpi)')
    ylabel('Molecules')
    hold off
    ax = gca;
    ax.FontSize = 20;
end

figure
yyaxis left
plot(Tv, P, '-k', 'Linewidth', 2)
set(gca, 'YColor', 'k')
ylabel('Degradation probability')

yyaxis right
plot(Tv, SiT(36,:), '-b', Tv, SiT(38,:), '--b', 'Linewidth', 2)
set(gca, 'YColor', 'b', 'YScale', 'log')
ylabel('Molecules')
xlabel('Time (HPI)')

legend('RNA degradation probability', 'Polymerase', 'NP protein', ...
       'FontSize', 13, 'Location', 'SE')

xticks([0 2 4 6 8 10 12])
ax = gca;
ax.FontSize = 20;


