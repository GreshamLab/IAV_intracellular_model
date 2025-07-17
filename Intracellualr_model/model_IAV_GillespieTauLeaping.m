function [SiT,Tv,SiTM,SiMC,Virs_G,Virs_C,Virs_time,NumS,step,NumG,mRNA,Prots,nRNP,pRNP,PM_cyt,ti,br] = model_IAV_GillespieTauLeaping(MOI,IniV,tf,printFreq,maxVir,maxVir_store,maxStep,maxCycle,e,Si,SiM,SiMC,ti,V0,Frac_I,Params,Reg,MR,MS,MP,NPs,dtMin)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script performs a stochastic simulation of Influenza A Virus (IAV) intracellular replication
% using the Gillespie algorithm with tau-leaping approximation. 

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

% Original stochastic model can be found in Heldt et al. 2015. We used a QSSA approximation to
% eliminate the free RNA and RNA bound to polymerases only, having only
% RNPs as states for the genomes and anti-genomes. That for, we used adapt
% the rate of the transcription of positive and negative genomes depending
% on the rate of the binding of RNA pol and N proteins. The same approach
% is used to eliminate the RNPs bound to M1 protein. 


% ============================
% ADD DEGRADATION REACTIONS
% ============================
% This section constructs and appends generic degradation reactions for each molecular species.
% These reactions are automatically generated and do not need to be included in the input matrices.
% However, the degradation rate constants for each species must be defined in Sheet 4 (parameter list).

% Get number of species from stoichiometric matrix
sM = size(MR);
NumS = sM(1,1);

% Create a vector indicating that all species undergo degradation
vDegs = ones(1, NumS);

% Build degradation stoichiometry:
% Each degradation reaction reduces one copy of its target species
MDegs1 = diag(vDegs) * -1;  % Net loss (stoichiometric matrix)
MDegs2 = diag(vDegs);       % Reactants (used in both MS and MP)

% Prepare new matrices with space for degradation reactions
sMR = size(MR);
sMP = size(MP);

% Expand stoichiometric matrix (MR) to include degradation reactions
MRwD = zeros(sMR(1), sMR(2) + NumS);
MRwD(:, 1:sMR(2)) = MR;                             % Original reactions
MRwD(:, sMR(2)+1:end) = MDegs1;                     % Degradation reactions

% Expand reactants matrix (MS) accordingly
MSwD = zeros(sMR(1), sMR(2) + NumS);
MSwD(:, 1:sMR(2)) = MS;
MSwD(:, sMR(2)+1:end) = MDegs2;

% Expand parameter matrix (MP) to include degradation rate constants
MPwD = zeros(sMP(1) + NumS, sMP(2) + NumS);
MPwD(1:sMP(1), 1:sMP(2)) = MP;
MPwD(sMP(1)+1:end, sMP(2)+1:end) = MDegs2;

% Replace original matrices with the expanded versions
MR = MRwD; 
MS = MSwD; 
MP = MPwD;

% Update matrix sizes
sM = size(MR);
NumS = sM(1,1);  % Total number of species (unchanged)
NumR = sM(1,2);  % Total number of reactions (including degradation)


% ======================
% MODEL DIMENSIONS AND INDICES
% ======================

% Number of protein types per virion (from NPs matrix)
NumP = size(NPs); 
NumP = NumP(1,2);  % Number of different proteins in the virion

% Number of IAV genomic segments (IAV has 8 RNA segments)
NumG = 8;

% Number of RNA species per segment (e.g., nRNA, pRNA, mRNA, and Matrix-bound cytoplasmic vRNA)
numRNAs = 4;

% Calculation of the virion degradation rate in the endosomes, based on the
% rate endosome espaces (Param(2)) and the fraction of successful ones (Params(3))
Params(3) = Params(2)*(1 - Params(3)) / Params(3);

% ======================
% STATE VECTOR COORDINATES
% ======================

% Each RNA species is organized in blocks per genome segment (NumG = 8),
% and each block contains 'numRNAs' variables.
% This computes the starting index of each segment block in the state vector.
V_coord = linspace(0, numRNAs * NumG, NumG + 1);  
V_coord = V_coord(1:NumG);  % Exclude the last endpoint

% Define the indices of different species types per segment:
nRNP     = V_coord + 4;  % Negative-sense ribonucleoproteins (nRNPs)
mRNA     = V_coord + 5;  % Messenger RNAs
pRNP     = V_coord + 6;  % Positive-sense replicative intermediate (pRNPs)
PM_cyt   = V_coord + 7;  % Matrix-bound cytoplasmic RNA species (for packing and export)

% ======================
% GLOBAL INDICES FOR KEY SPECIES
% ======================

% Indices of viral proteins in the state vector
Prots = 36:42;  % General range for proteins

% Specific proteins of interest:
Pol = 36;  % RNA polymerase complex
NP  = 38;  % Nucleoprotein
M1  = 40;  % Matrix protein 1
NEP = 42;  % Nuclear export protein

% Other index groups for grouping and processes:
RNA_nuc       = 4;         % Nuclear importation
Pos_syn       = 5:12;      % Positive-strand synthesis (per segment)
Transcriptions = 13:20;    % mRNA transcription events (per segment)
Neg_syn       = 21:28;     % Negative-strand synthesis (per segment)
MBind         = 29:36;     % Matrix protein binding events (e.g., for export/packaging)
Translations  = 37:43;     % Protein translation events (per gene product)
Vir_Prod      = 44;        % Virion production event index

% Grouped indices for synthesized RNAs
Gen_syn = [Pos_syn; Neg_syn];  % Group both strands' synthesis events for convenience

% ============================
% REACTION INITIALIZATION
% ============================

% Preallocate matrices for reaction parameters and bookkeeping

MParsM = zeros(MOI, NumR);  % Matrix to store parameter values for each reaction and phenotype (MOI variants)
MPars  = zeros(1, NumR);    % Final effective parameter values per reaction (averaged over MOI variants)
RegM   = zeros(MOI, length(Reg));     % Parameters of reactions excluded via QSSA, one row per variant
IndS   = zeros(3, NumR);    % Matrix to store indices of up to 3 reactants per reaction
gi     = ones(1, NumS);     % Vector to store the maximum order (number of reactants) per species
Rorder = ones(1, NumR);     % Vector with the number of reactants per reaction

% ============================
% Identify reactants per reaction
% ============================

for J = 1:NumR
    I = find(MS(:, J));         % Find indices of species that are reactants in reaction J
    s = size(I); s = s(1);      % Number of reactants in this reaction
    IndS(1:s, J) = I;           % Store reactant indices
    Rorder(1, J) = s;           % Store reaction order (number of reactants)
end

% Replace empty entries in IndS with NumS+1 (a dummy index outside the real species list)
IndS(IndS == 0) = NumS + 1;

% Set order to 1 for reactions with no reactants (prevents division by zero or indexing errors)
Rorder(Rorder == 0) = 1;

% Set translation reactions (first of their group) to order 1 (could relate to special handling)
Rorder(Translations(1)) = 1;

% ============================
% Compute max reaction order per species
% ============================

for S = 1:NumS
    [aI, bI] = find(IndS == S);  % Find all reactions where species S is a reactant
    gi(S) = max(Rorder(bI));    % Record the highest order in which this species participates
end

% Default any missing values to 1 (for robustness)
gi(gi == 0) = 1;

% ============================
% INITIAL CONDITIONS PER VARIANT (MOI)
% ============================

for m = 1:MOI
    % Set initial amount of virions for species 1 (e.g., initial infection input)
    Si(1,1)     = MOI * IniV;  % Total initial virions
    SiM(m,1)    = IniV;        % Virions of variant m
    SiMC(m,1)   = IniV;        % Same value, used across cycles (if applicable)

    % Scale MP by parameter values; here, noise could be added (commented out)
    MP2 = MP .* Params; % .* (0.5 + rand(sP,1));  % Element-wise scaling of kinetic matrix

    % Copy QSSA-eliminated parameters to RegM for this variant
    RegM(m,:) = Reg; % .* (0.5 + rand(4,1));  % Option to add stochastic variation
    
    % Compute the effective rate constant for each reaction:
    for J = 1:NumR
        idx = find(MP2(:, J));                         % Find parameter indices used in reaction J
        MParsM(m, J) = prod(MP2(idx, J));              % Multiply them to compute the total rate
    end
end

% ============================
% COMPUTE AVERAGE REACTION RATES OVER MOI VARIANTS
% ============================

% These are dummy initial species vectors, normalized across MOI variants
Si1  = ones(1, NumS + 1);
SiM1 = ones(MOI, NumS + 1) ./ MOI;

for J = 1:NumR
    % mS contains normalized concentration of the limiting reactant across MOI variants
    mS = SiM1(:, IndS(Rorder(J), J));  % Get amount of highest-order reactant in reaction J
    mS = mS ./ sum(mS);                % Normalize across variants
    mS(isnan(mS)) = 1 / MOI;           % Handle zero-division by assigning uniform probability

    % Weighted average of the rate constant across variants (based on mS weights)
    MPars(1, J) = sum(mS .* MParsM(:, J));
end
 
% ============================
% PROPENSITY VECTORS AND OUTPUT STORAGE
% ============================

% Number of Gillespie steps to simulate
nSteps = 1e5;

% Initialize time vector for each simulation step
Tv = zeros(1, nSteps);  % Time points

% Initialize matrices to store species counts at each time point
SiT  = zeros(NumS, nSteps);         % State vector for global (aggregated) simulation
SiTM = zeros(NumS * MOI, nSteps);   % State vectors per MOI variant (stacked)

% Initialize tracking of genotypes present in virions produced
Frac_G    = zeros(NumG, maxVir * 2);  % Stores the initial virion of each segment packed
Frac_C    = zeros(NumG, maxVir * 2);  % Stores the number of replication cycles of each segment packed
Vir_time  = zeros(1, maxVir * 2);     % Stores time of virion production
nv = 0; % Number of stored virions

% ============================
% REACTION INDEXING AND MASKING
% ============================

% Create index matrix mapping each cytoplasmic PM RNA entry (from SiM) to its variant
IndGen = repmat(1:size(SiM(:, PM_cyt), 1), 1, size(SiM(:, PM_cyt), 2));
IndGen = reshape(IndGen, size(SiM(:, PM_cyt)));

% Initialize simulation step counter
step = 0;

% MRM: Modified reaction matrix used for computing net change in subset of reactions
MRM = MR;

% Mask out selected species contributions in specific reactions
MRM(M1, MBind(:))         = 0;  % Ignore M1 contribution in matrix binding reactions
MRM(NEP, MBind(:))        = 0;  % Ignore NEP contribution in matrix binding reactions
MRM([Pol NP], Gen_syn(:)) = 0;  % Ignore polymerase and NP contributions in RNA synthesis
MRM(Prots, Vir_Prod)      = 0;  % Ignore protein contributions in virion production


tic
br = 0;
ap = 0;

%% ============================
%% SIMULATION LOOP
%% ============================

while (ti < tf && step < maxStep && nv < maxVir)

    step = step + 1;  % Advance simulation step
    
    % Copy current global and per-variant state vectors
    Si1(1:NumS) = Si;
    SiM1(:,1:NumS) = SiM;
    
    % Compute effective rate constants for each reaction across variants
    for J = 1:NumR
        mS = SiM1(:, IndS(Rorder(J), J));         % Abundance of key reactant in each variant
        mS = mS ./ sum(mS);                       % Normalize
        mS(isnan(mS)) = 1 / MOI;                  % Handle division by zero
        MPars(1, J) = sum(mS .* MParsM(:, J));    % Weighted average over variants
    end

    % Update QSSA-regulated parameters using weighted contributions from each variant

    mR = SiM1(:, Pol) ./ Si1(Pol); mR(isnan(mR)) = 1 / MOI;
    Reg(3) = sum(RegM(:,3) .* mR);  % Polymerase-regulated

    mR = SiM1(:, NP) ./ Si1(NP); mR(isnan(mR)) = 1 / MOI;
    Reg(4) = sum(RegM(:,4) .* mR);  % NP-regulated

    mR = SiM1(:, NEP) ./ Si1(NEP); mR(isnan(mR)) = 1 / MOI;
    Reg(6) = sum(RegM(:,6) .* mR);  % NEP-regulated

    % ----------------------------
    % Compute reaction propensities
    % ----------------------------
    A1 = Propensities(Si1, IndS, MPars, Vir_Prod, Params, PM_cyt, Prots, NPs, ...
                             Translations, mRNA, Transcriptions, Pol, Neg_syn, ...
                             Pos_syn, Reg, NP, MR, Si, MBind, NEP, SiM1, NumS);

    % Compute next time step and critical reactions
    [dt, dt_noC, dt_C, Rcrit, Rnocrit] = dtCal(Si, PM_cyt, MR, Vir_Prod, ...
                                                     NumR, MS, A1, MR, e, gi);

    if dt > tf  % If next step exceeds final time, exit
        step = step - 1; 
        break; 
    end

    % Approximate fast reactions if timescales collapse and we're past early transient
    if min([dt_C dt_noC]) < dtMin && dt_C > 0 && dt_noC > 0 && ti > 2*60
        ap = ap + 1;
        [Si, SiM, SiMC, A1, dt, dt_noC, dt_C, Rcrit, Rnocrit] = ...
            nRNPapprox(Si, SiM, SiMC, A1, MR, SiM1, nRNP, PM_cyt, Params, ...
                             NumS, maxCycle, IndS, MBind, Neg_syn, Prots, M1, Pol, ...
                             NP, NEP, Reg, Vir_Prod, NumR, MS, e, gi, MOI, ...
                             MParsM, NP, dt, dt_noC, dt_C);
    end

    % Sanity checks
    if sum(Si < 0) > 0
        Si, ti, SiM  % Debug output if negative values occur
    end
    if sum(A1) == 0
        step = step - 1; 
        break; 
    end

    % ----------------------------
    % Gillespie or Tau-Leap execution
    % ----------------------------
    if dt_noC < (dt * 10)  % Condition to use full Gillespie
        [Si, SiM, SiMC, Frac_G, Frac_C, nv, Vir_time] = ...
            model_gillespie(Si, SiM, SiMC, Si1, SiM1, A1, dt, MR, MRM, IndS, ...
                            Pos_syn, RNA_nuc, V0, nRNP, 1:NumR, Frac_G, Frac_C, ...
                            nv, Vir_time, Vir_Prod, NPs, PM_cyt, Prots, NumG, ...
                            NumS, maxCycle, ti, Frac_I, MBind, M1, MParsM, Pol, ...
                            NP, NEP, Reg, NumP, IndGen);
    else
        if dt_noC < dtMin
            dt_noC = dtMin;  % Enforce minimum step size
        end
        dt = min([dt_noC, dt_C]);  % Choose conservative step size

        % Tau-leap step (non-critical reactions)
        [Si, SiM, SiMC, Frac_G, Frac_C, nv, Vir_time] = ...
            model_tauL(Si, SiM, SiMC, Si1, SiM1, A1, dt, MR, MRM, IndS, ...
                       Pos_syn, RNA_nuc, V0, nRNP, Rnocrit, Frac_G, Frac_C, ...
                       nv, Vir_time, Vir_Prod, NPs, PM_cyt, Prots, NumG, ...
                       NumS, maxCycle, ti, Frac_I, MBind, M1, MParsM, Pol, ...
                       NP, NEP, Reg, NumP, IndGen);

        % Gillespie for critical reactions if needed
        if dt_C <= dt
            [Si, SiM, SiMC, Frac_G, Frac_C, nv, Vir_time] = ...
                model_gillespie(Si, SiM, SiMC, Si1, SiM1, A1, dt, MR, MRM, IndS, ...
                                Pos_syn, RNA_nuc, V0, nRNP, Rcrit, Frac_G, Frac_C, ...
                                nv, Vir_time, Vir_Prod, NPs, PM_cyt, Prots, NumG, ...
                                NumS, maxCycle, ti, Frac_I, MBind, M1, MParsM, ...
                                Pol, NP, NEP, Reg, NumP, IndGen);
        end
    end

    % Advance simulation time
    ti = ti + dt;

    % Print progress every printFreq steps
    if mod(step, printFreq) == 0
        ti / 60      % Current time (in hours)
        sum(Si(NumS)) % Number of complete virions
        toc
        tic
    end

    % Consistency check: verify that the global state equals the sum over variants
    if prod(Si(1:end-1) == sum(SiM(:, 1:end-1), 1)) == 0
        step = step - 1; 
        break; 
    end

    % Store output
    SiT(:, step) = Si;          % Global species counts
    Tv(1, step) = ti / 60;      % Time vector (in hours)
    for m = 1:MOI
        IniS = (m - 1) * NumS + 1;
        FinS = m * NumS;
        SiTM(IniS:FinS, step) = SiM(m, :);  % Per-variant states
    end

    % Emergency stop: 2 hours of runtime exceeded
    A = toc;
    if A > 60 * 60 * 2
        br = 1;
        step = step - 1;
        break;
    end

    % Emergency stop: system essentially dead (no genomes, no replication)
    if any((Si(nRNP) + Si(pRNP) + Si(PM_cyt)) < 0.1) && sum(Si(1:3)) < 0.1
        Si  % Final state output
        break;
    end
end

%%
% If the recorded time points exceed the preallocated maximum (`step`), truncate the arrays
% This prevents oversized output due to unexpected extra steps
if length(Tv) > step
    SiT = SiT(:,1:step);      % Truncate total species time trace
    Tv = Tv(1:step);          % Truncate time vector
    SiTM = SiTM(:,1:step);    % Truncate multi-variant species trace
end

% If the number of virions produced (`nv`) exceeds the max allowed for output (`maxVir_store`)
% randomly subsample a representative set of virions
if nv > maxVir_store
    random_vir = randsample(nv, maxVir_store);     % Randomly sample indices
    Virs_G = Frac_G(:, random_vir);                % Sampled genome identities
    Virs_C = Frac_C(:, random_vir);                % Sampled cycle identifiers
    Virs_time = Vir_time(random_vir);              % Corresponding production times
else
    % If nv is below threshold, pad output by selecting the first `maxVir_store` entries
    Virs_G = Frac_G(:, 1:maxVir_store);
    Virs_C = Frac_C(:, 1:maxVir_store);
    Virs_time = Vir_time(1:maxVir_store);
end


end

