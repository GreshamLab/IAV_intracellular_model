%% IAV STOCHASTIC SIMULATION – MOI SWEEP
% -------------------------------------------------------------------------
% This script runs stochastic simulations of Influenza A Virus (IAV)
% intracellular replication for a range of MOIs (Multiplicity of Infection).
% For each MOI, 1000 replicate simulations are performed.
%
% Required input matrices are loaded from an Excel file.
% For each replicate, the script initializes state variables, runs the
% simulation using a Gillespie/tau-leaping solver, and stores output data.
%
% Results are saved as CSV files in a time-stamped output directory.
% -------------------------------------------------------------------------

% ===============================================================================
% This script is used to generate simulated data like the ones used 
% for the Figure 2A and Figure 3A-B of the paper Segredo-Otero and Gresham 2025.
% ===============================================================================

clear all
close all
clc

%% CREATE OUTPUT DIRECTORY ------------------------------------------------
path = pwd;  % Get current working directory

% Generate timestamped folder name for output
date = sprintf(datestr(now,'mm-dd-yyyy_HH-MM-SS'));
dirName = strcat('Result_', date);
mkdir(path, dirName);

Path = fullfile(path, dirName);  % Full output path

%% LOAD INPUT EXCEL WORKBOOK ----------------------------------------------
% The workbook **Param_mat** must contain at least six sheets:
%   1. Stoichiometric matrix  (species × reactions)
%   2. Reactants matrix       (species × reactions, logical)
%   3. Rate‑constant matrix   (parameters × reactions)
%   4. Parameter list         (symbols + values)
%   5. QSSA‑eliminated params (fast reactions handled analytically)
%   6. Virion composition     (protein copies per virion)
Param_mat = 'Stock_matrix_MDCK.xlsx';

MR     = xlsread(Param_mat, 1);      % Stoichiometric matrix
MS     = xlsread(Param_mat, 2);      % Reactants matrix
MP     = xlsread(Param_mat, 3);      % Rate constants
Params = xlsread(Param_mat, 4);      % Parameters
Reg    = xlsread(Param_mat, 5);      % QSSA-regulated parameters
NPs    = xlsread(Param_mat, 6)';     % Virion structure

%% TIME STEP SETTINGS FOR EACH MOI ----------------------------------------
dtM = zeros(1,20);           % Preallocate vector of dt values

% Use smaller dt for low MOI (to increase stability)
dtM(1,1:8)  = 0.001;
dtM(1,9:20) = 0.02;

%% LOOP OVER MOI ----------------------------------------------------------
for MOI = 1:20

    % Parameters constant across replicates
    Rep = 1000;           % Number of simulation replicates
    printFreq = 5e5;      % Console printing frequency
    dtLim = 0.2;          % Max allowed timestep
    maxVir = 1e5;         % Max virions per replicate
    maxVir_store = 1e3;   % Max virions to retain per replicate
    maxStep = 1e5;        % Max Gillespie steps
    maxCycle = 30;        % Max number of replication cycles stored
    e = 0.03;             % Tau-leap error tolerance
    tf = 10 * 60;         % Simulation time (10 hours)
    dtMin = dtM(MOI);     % Minimum timestep
    IniV = 1;             % Virions per genotype
    NumS = size(MR,1);    % Number of species
    NumG = 8;             % Number of genome segments

    % Preallocate outputs
    Fit = zeros(1, Rep);
    Tit_int = 1000;                       % Output timepoints to store
    Tit = zeros(Tit_int, Rep);           % Timecourse of virion counts
    Times = zeros(Tit_int, Rep);         % Time vector for Tit
    Time = zeros(1, Rep);                % Time at end of sim
    GENS = zeros(NumG, MOI, Rep);        % Genome count per segment
    Virs_G = zeros(NumG, maxVir_store, Rep);
    Virs_C = zeros(NumG, maxVir_store, Rep);
    Virs_time = zeros(maxVir_store, Rep);
    Muts_per_virion = zeros(1, Rep);     % Placeholder for future use

    % Show current MOI
    disp(['Running MOI = ', num2str(MOI)]);

    for rep = 1:Rep

        % Initialize state vectors
        Si   = zeros(1, NumS);                  % Global state
        SiM  = zeros(MOI, NumS);                % Per-genotype state
        SiMC = zeros(MOI, NumS * maxCycle);     % With replication cycles
        ti   = 0;                               % Initial time

        % Initialize virion genome composition (complete segments)
        V0 = ones(MOI, NumG);

        % Initialize fraction of infecting virions by genotype class
        Frac_I = zeros(1,8); 
        Frac_I(1) = 0.8;   % 80% complete
        Frac_I(2) = 0.2;   % 20% missing one segment

        % Run simulation
        [SiT, Tv, SiTM, SiMC, Frac_G, Frac_C, Virions_time, ...
             NumS, step, NumG, mRNA, Prots, nRNP, pRNP, PM_cyt, ...
             ti, br] = model_IAV_GillespieTauLeaping(MOI, IniV, tf, ...
             printFreq, maxVir, maxVir_store, maxStep, ...
             maxCycle, e, Si, SiM, SiMC, ti, V0, Frac_I, ...
             Params, Reg, MR, MS, MP, NPs, dtMin);

        % Reshape time series to (species × time × genotype)
        SiTM1 = zeros(NumS,step,MOI);
        for m = 1:MOI
            IniS = (m-1)*NumS + 1;
            FinS = m*NumS;
            SiTM1(:,:,m) = SiTM(IniS:FinS,:);
        end
        SiTM = SiTM1;

        % Retrieve final genome counts per genotype
        SiM = reshape(SiTM(:,end,:),NumS,MOI); SiM = SiM';
        for m = 1:MOI
            GENS(:,m,rep) = floor(SiM(m,PM_cyt));
        end

        % Final number of virions produced
        Fit(1,rep) = SiT(NumS,end);
        Time(1,rep) = Tv(step);

        % Save evenly spaced timepoints from simulation
        Idx = round(linspace(1,step,Tit_int));
        Tit(:,rep) = SiT(NumS,Idx);
        Times(:,rep) = Tv(Idx);

        % Save virion composition info
        Virs_G(:,:,rep) = Frac_G;
        Virs_C(:,:,rep) = Frac_C;
        Virs_time(:,rep) = Virions_time;
    end

    %% EXPORT TO CSV ------------------------------------------------------

    % Reshape outputs
    Gens1 = reshape(GENS,NumG*Rep,MOI);
    Virs_G1 = reshape(Virs_G,NumG*maxVir_store,Rep);
    Virs_C1 = reshape(Virs_C,NumG*maxVir_store,Rep);
    Fit1 = [Fit;Time];

    % Write to CSV files
    moiT = num2str(MOI*IniV);
    base = strcat('MOI_',moiT,'_');

    writematrix(Gens1,      fullfile(Path, [base, 'Genomes_', date, '.csv']));
    writematrix(Fit1,       fullfile(Path, [base, 'Fit_', date, '.csv']));
    writematrix(Virs_G1,    fullfile(Path, [base, 'Virs_G_', date, '.csv']));
    writematrix(Virs_C1,    fullfile(Path, [base, 'Virs_C_', date, '.csv']));
    writematrix(Virs_time,  fullfile(Path, [base, 'Virs_time_', date, '.csv']));
    writematrix(Tit,        fullfile(Path, [base, 'Tit_', date, '.csv']));
    writematrix(Times,      fullfile(Path, [base, 'Times_', date, '.csv']));
    
end

%% CLOSE PARALLEL POOL IF ACTIVE ------------------------------------------
parpool().delete;
