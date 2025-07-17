%% POPULATION‑LEVEL IAV INFECTION SIMULATION COMPARISON WITH EXPERIMENTAL DATA
% -------------------------------------------------------------------------
% This script uses pre‑generated single‑cell infection results to build a
% population‑scale simulation (in this case, we're using
% the ones form Result_03-23-2025_08-14-13, but these can be generated with 
% the script Central_IAVGillespie_multiple). It compares a neutral "null model" virus (all
% infections succeed, identical output per virion infecting cells) with 
% Influenza A Virus (IAV) simulation results. The simulation tracks 2048 neutral
% barcode variants, as donde in  Varble et al. 2014., computes their 
% frequencies before/after infection, and estimates effective population size (Ne).
% The experimental data form Varble et al. 2014 used by this script is
% provided in Data.xlsx
% -------------------------------------------------------------------------

% =====================================================================
% This script is used to generate simulated data like the ones used 
% for the Figure 4B-c of the paper Segredo-Otero and Gresham 2025.
% =====================================================================


clear all
close all
clc;

%% LOAD PRE‑COMPUTED SINGLE‑CELL RESULTS ----------------------------------
current = strcat(pwd, '/Result_03-23-2025_08-14-13');  % Folder with CSVs
path    = current(end-19:end);   % Timestamp substring used in filenames

M      = 20;      % Max MOI simulated
MOIs   = 1:M;     % Vector 1…20
NumG   = 8;       % Genome segments
Rep    = 1000;    % Replicates per MOI
maxVir = 1e3;     % Virions stored per replicate

% Preallocate storage matrices
Genomes_store = nan(NumG*Rep, M, M);
Fit_store     = nan(Rep, M);

for m = 1:M
    % Construct filenames for current MOI
    Genomes_name = strcat(current, '/MOI_', num2str(MOIs(m)), '_Genomes', path, '.csv');
    Fit_name     = strcat(current, '/MOI_', num2str(MOIs(m)), '_Fit',     path, '.csv');

    % Read CSVs
    Genomes = csvread(Genomes_name);
    Fit     = csvread(Fit_name);

    % Store
    Genomes_store(:, 1:MOIs(m), MOIs(m)) = Genomes;
    Fit_store(:, MOIs(m)) = Fit(1, :);
end

%% PROCESS GENOME & FITNESS DATA PER MOI ---------------------------------
mFreqAll = zeros(M, M);  % Mean genome frequencies per segment per MOI
mFit     = zeros(1, M);  % Mean virion output per replicate
vFreqAll = zeros(M, M);  % SD of frequencies
vFit     = zeros(1, M);  % SD of virion output
mFreq    = zeros(1, M);  % Fraction of successful infections (>10 virions)
vFreq    = zeros(1, M);

for m = 1:M
    Genomes1 = reshape(Genomes_store(:,1:m,m), NumG, m, Rep);
    Fit      = Fit_store(:, m);

    FreqAll = zeros(M, Rep);  % Segment‑wise frequencies per replicate

    for rep = 1:Rep
        Gens = Genomes1(:,:,rep);
        if Fit(rep) > 0  % Replicate produced virions
            A = mean(sort(Gens ./ sum(Gens,2), 2, 'descend'));  % Mean segment freq
            FreqAll(1:m, rep) = A;
        end
    end

    % Success fraction & variance
    mFreq(m) = mean(Fit > 10);
    vFreq(m) = sqrt(var(Fit > 10));

    % Remove zero‑output replicates
    FreqAll(:, Fit == 0) = [];
    Fit(Fit == 0) = [];

    % Mean/SD of virion output & segment frequencies
    mFit(m) = mean(Fit);
    vFit(m) = sqrt(var(Fit));
    mFreqAll(:, m) = nanmean(FreqAll');
    vFreqAll(:, m) = sqrt(nanvar(FreqAll'));
end

%% LOAD POPULATION DATA FROM Varble et al. 2014 ---------------------------
Data_2014 = 'Data.xlsx';
ndata = xlsread(Data_2014, 1);
Total_differentes_seqs = length(ndata);
Ini_title = 1e4;              % Initial PFU in P0
Ini_freq  = 1 / Total_differentes_seqs;

%% POPULATION‑LEVEL FORWARD SIMULATION ------------------------------------
Reps = 100;                 % Population replicates
NE_vir = zeros(1, Reps);    % Ne for null model
NE_IAV = zeros(1, Reps);    % Ne for IAV model

tH = 3.2;                   % Virion half‑life [h] from Perelson et al. 2006
dV = -log(0.5) / tH;        % Decay rate

dt = 12;                    % Time per viral generation [h]
Titer_per_generation = [0 3e4 1e7 2e7 1e7] * 2; % Published PFU/plate in Varble et al. 2014

% Compute virion decay and production between generations
DegV    = abs(Titer_per_generation(1:end-1) .* (exp(-dV*dt) - 1));
dTiter  = diff(Titer_per_generation .* (1 + exp(-dV*dt/2)));

% Average virion output per cell (from single‑cell model)
Fit1  = Fit_store(:,1);  Fit20 = Fit_store(:,M);
mFit1 = mean(Fit1(Fit1 > 0));  mFit20 = mean(Fit20(Fit20 > 0));
Av_fit_per_cells = [mFit1 mFit20 mFit20 mFit20];

Frac_PFU = mFreq(1) * 0.8;   % PFU fraction (successful infections)
Cells_per_generation = round(((dTiter + DegV) ./ Frac_PFU) ./ Av_fit_per_cells);
MOI_vec = [1 M M M];        % MOI per generation
nG = 2 * 24 / 12;           % 2 generations/day → nG = 4

for reps = 1:Reps
    % Initial barcode frequencies
    Seqs_freq = ones(1, Total_differentes_seqs) * Ini_freq;
    Num_seqs  = mnrnd(Ini_title, Seqs_freq);
    Seqs_freq_vir = Num_seqs / sum(Num_seqs); % Null model
    Seqs_freq_IAV = Seqs_freq_vir;            % IAV model (same start)

    for n = 1:nG
        Cells = Cells_per_generation(n);
        cMOI  = MOI_vec(n);

        % Draw single‑cell infection outputs for current MOI
        Genomes   = reshape(Genomes_store(:,1:cMOI,cMOI), [], cMOI, Rep);
        Infections_MOI = Genomes(6,:,:)./sum(Genomes(6,:,:),2);
        Infections_MOI(isnan(Infections_MOI)) = 0;

        % Sample which barcodes infect each cell
        Cells_population_IAV = randsample(Total_differentes_seqs, cMOI*Cells, true, Seqs_freq_IAV);
        Cells_population_IAV = reshape(Cells_population_IAV, cMOI, Cells);
        Cells_population_vir = randsample(Total_differentes_seqs, cMOI*Cells, true, Seqs_freq_vir);
        Cells_population_vir = reshape(Cells_population_vir, cMOI, Cells);

        Rand_infs = randi(Rep, 1, Cells);  % Pick random replicate traces
        Infections_IAV = reshape(Infections_MOI(:,:,Rand_infs), cMOI, []) .* ...
                         Fit_store(Rand_infs, cMOI)';

        Seqs_count_vir = zeros(1, Total_differentes_seqs);
        Seqs_count_IAV = zeros(1, Total_differentes_seqs);

        for m = 1:cMOI
            % Null model: each infection contributes 1 virion
            Seqs_count_vir = Seqs_count_vir + histcounts(Cells_population_vir(m,:), 1:Total_differentes_seqs+1);
            % IAV model: contribution weighted by single‑cell output
            Seqs_count_IAV = Seqs_count_IAV + ...
                accumarray(Cells_population_IAV(m,:)', Infections_IAV(m,:)', ...
                           [Total_differentes_seqs, 1], @sum, 0)';
        end

        % Update frequencies
        Seqs_freq_vir = Seqs_count_vir / sum(Seqs_count_vir);
        Seqs_freq_IAV = Seqs_count_IAV / sum(Seqs_count_IAV);
    end

    % Effective population size estimates
    NE_vir(reps) = mean(Seqs_freq .* (1-Seqs_freq)) / mean((Seqs_freq_vir - Seqs_freq).^2);
    NE_IAV(reps) = mean(Seqs_freq .* (1-Seqs_freq)) / mean((Seqs_freq_IAV - Seqs_freq).^2);
end

%% COMPUTE Ne FROM VARBLE ET AL. 2014 DATA --------------------------------
Seqs_freq_data      = ndata(:,3);
Seqs_freq_IAV_data  = ndata(:,4:6);
Ne_IAV_data = mean(Seqs_freq_data .* (1-Seqs_freq_data)) ./ ...
              mean((Seqs_freq_IAV_data - Seqs_freq_data).^2);

% Means and SDs
mNe_data = mean(Ne_IAV_data);
vNe_data = sqrt(var(Ne_IAV_data,0,2));
mNe_vir  = mean(NE_vir);
vNe_vir  = sqrt(var(NE_vir,0,2));
mNe_IAV  = mean(NE_IAV);
vNe_IAV  = sqrt(var(NE_IAV,0,2));

mNEs = [mNe_data mNe_IAV mNe_vir];
vNEs = [vNe_data vNe_IAV vNe_vir];

%% PLOT EFFECTIVE POPULATION SIZE ----------------------------------------
figure('Name','Effective Population Size'); hold on
bar(mNEs)
errorbar(1:3, mNEs, vNEs, '*k')
xticks(1:3)
xticklabels({'Data','IAV model','Null model'})
set(gca,'YScale','log','YLim',[1e1 3e3])
ax = gca; ax.FontSize = 20;

%% PLOT – FREQUENCY SHIFTS BEFORE VS AFTER INFECTION -------------------------
figure('Name','Frequency Shift – Example'); hold on
antes   = Seqs_freq;       % Initial barcode frequencies
despues = Seqs_freq_IAV;   % Final barcode frequencies after IAV model

numMutantes = length(antes);      % Number of tracked variants
colors = rand(numMutantes, 3);    % Assign random color to each variant

% Bar plot settings
barWidth  = 1;
xAntes    = 1;
xDespues  = 3;

% Plot stacked bars for initial (P0) and final (P1) frequencies
hAntes = bar(xAntes, antes', 'stacked', 'EdgeColor', 'none', 'BarWidth', barWidth);
hDespues = bar(xDespues, despues', 'stacked', 'EdgeColor', 'none', 'BarWidth', barWidth);

% Assign consistent color per variant across bars
for i = 1:numMutantes
    hAntes(i).FaceColor = colors(i, :);
    hDespues(i).FaceColor = colors(i, :);
end

% Compute cumulative heights of stacked segments
yAntes   = [0, cumsum(antes)];
yDespues = [0, cumsum(despues)];

% Draw connecting polygons between segments ("alluvial flow")
for i = 1:numMutantes
    xCoords = [xAntes+barWidth/2, xDespues-barWidth/2, xDespues-barWidth/2, xAntes+barWidth/2];
    yCoords = [yAntes(i), yDespues(i), yDespues(i+1), yAntes(i+1)];
    fill(xCoords, yCoords, colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end

xlim([xAntes-barWidth/2, xDespues+barWidth/2]);
ylim([0, 1]);
xticks([xAntes, xDespues]);
xticklabels({'P0', 'P1'});
ylabel('Frequencies');
ax = gca; ax.FontSize = 15;
hold off;
