%% FITNESS–MOI PARAMETER SWEEP IN IAV POPULATION SIMULATION
% -------------------------------------------------------------------------
% This script performs a grid simulation sweeping across a range of fitness 
% advantages and multiplicities of infection (MOI). It uses single-cell infection 
% results from IAV simulations to model population-level dynamics (in this case, we're using
% the ones form Result_03-23-2025_08-14-13, but these can be generated with 
% the script Central_IAVGillespie_multiple).
% For each parameter combination, the simulation introduces a single high-fitness 
% variant and tracks its final frequency after 14 days (~28 generations). 
% Results are summarized in a heatmap showing how selection and MOI affect 
% variant fixation probability.
% -------------------------------------------------------------------------

% =====================================================================
% This script is used to generate simulated data like the ones used 
% for the Figure 5B of the paper Segredo-Otero and Gresham 2025.
% =====================================================================

clear all
close all
clc

%% Load precomputed single-cell infection results ---------------------------

current = strcat(pwd,'/Result_03-23-2025_08-14-13');
path = current(end-19:end);
M = 20;                  % Maximum MOI value in dataset
MOIs = 1:M;
NumG = 8;                % Number of genomes per simulation
Rep = 1000;              % Number of replicates
Genomes_store = nan(NumG*Rep,M,M); % Storage for genomes per MOI and replicate
Fit_store = nan(Rep,M);             % Storage for fitness values per MOI

for m = 1:M
    % Load genomes and fitness data for each MOI
    Genomes_name = strcat(current,'/MOI_',num2str(MOIs(m)),'_Genomes',path,'.csv');
    Fit_name = strcat(current,'/MOI_',num2str(MOIs(m)),'_Fit',path,'.csv');
    Genomes = csvread(Genomes_name);
    Fit = csvread(Fit_name);
    
    % Clean and store
    Genomes_store(:,1:MOIs(m),MOIs(m)) = Genomes;
    Fit = Fit(1,:);
    Fit(isnan(Fit)) = 0; Fit(Fit == Inf) = 0;
    Fit_store(:,MOIs(m)) = Fit;
end

%% Summarize frequency and fitness across MOIs ------------------------------

mFreqAll = zeros(M,M); vFreqAll = zeros(M,M);
mFit = zeros(1,M); vFit = zeros(1,M);
mFreq = zeros(1,M); vFreq = zeros(1,M);

for m = 1:M
    Genomes1 = reshape(Genomes_store(:,1:m,m),NumG,m,Rep);
    Fit = Fit_store(:,m);
    FreqAll = zeros(M,Rep);

    for rep = 1:Rep
        Gens = Genomes1(:,:,rep);
        if Fit(rep) > 0
            A = mean(sort(Gens./sum(Gens,2),2,'descend'));
            if m < M
                FreqAll(1:m,rep) = A;
            else
                FreqAll(:,rep) = A;
            end
        end
    end

    mFreq(1,m) = mean(Fit>10); vFreq(1,m) = sqrt(var(Fit>10));
    FreqAll(:,Fit==0) = [];
    Fit(Fit==0) = [];
    mFit(1,m) = mean(Fit); vFit(1,m) = sqrt(var(Fit));
    mFreqAll(:,m) = nanmean(FreqAll'); vFreqAll(:,m) = sqrt(nanvar(FreqAll'));
end

%% Set simulation parameters ------------------------------------------------

Total_differentes_seqs = 1e4;        % Number of barcode variants
Ini_title = 1e4;                     % Initial total barcode counts
Ini_freq = 1/Total_differentes_seqs; % Initial frequency per variant
Reps = 100;                          % Repeats per condition

Rep_fit = 10;                        % Number of fitness values tested
W_factor = linspace(1,5,Rep_fit);    % Fitness of mutant (w > 1)
vMOI = linspace(1,10,Rep_fit);       % MOI values tested

Days = 14;                           % Number of days
nG = Days*2;                         % 2 generations per day

F1_IAV = zeros(Rep_fit,Rep_fit,Reps);    % Stores final freq of mutant
F1_IAV(1,:,:) = Ini_freq;               % Initial freq = 1/N
Segment = randi(8);                     % Random segment from 1–8 for infection profiles

%% Run grid simulation: sweep over MOI and fitness --------------------------

for rep_MOI = 1:Rep_fit
for rep_fit = 1:Rep_fit

    MOI = vMOI(rep_MOI);
    w_adv = W_factor(rep_fit);
    Cells_per_generation = ones(1,nG).*round(1e5/MOI); % Keep total virions constant

    w_factor = ones(1,Total_differentes_seqs); w_factor(1) = w_adv; % Mutant at index 1
    w_factor = w_factor./sum(w_factor); % Normalize to keep total output constant

    for reps = 1:Reps
        % Reset initial barcode distribution
        Seqs_freq = ones(1,Total_differentes_seqs).*Ini_freq;
        Num_seqs = mnrnd(Ini_title,Seqs_freq);
        Seqs_freq_IAV = Seqs_freq;

        % Simulate nG generations
        for n = 1:nG
            Cells = Cells_per_generation(n);
            cMOI = MOI;
            Genomes = reshape(Genomes_store(:,1:cMOI,cMOI),[],cMOI,Rep);
            Infections_MOI = Genomes(Segment,:,:)./sum(Genomes(Segment,:,:),2);
            Infections_MOI(isnan(Infections_MOI)) = 0;

            % Simulate cell population
            Cells_population_IAV = randsample(Total_differentes_seqs,cMOI*Cells,true,Seqs_freq_IAV);
            Cells_population_IAV = reshape(Cells_population_IAV,cMOI,Cells);
            Rand_infs = randi(Rep,1,Cells);
            Infections_IAV = reshape(Infections_MOI(:,:,Rand_infs),cMOI,[]).*Fit_store(Rand_infs,cMOI)';

            Seqs_count_IAV = zeros(1,Total_differentes_seqs);
            for m = 1:cMOI
                Seqs_count_IAV(1,:) = Seqs_count_IAV(1,:) + ...
                    (accumarray(Cells_population_IAV(m,:)', Infections_IAV(m,:)', ...
                    [Total_differentes_seqs, 1], @sum, 0))' .* w_factor;
            end

            % Normalize to get new frequency
            Seqs_freq_IAV = Seqs_count_IAV./sum(Seqs_count_IAV);
        end

        % Store final frequency of the beneficial mutant (index 1)
        F1_IAV(rep_MOI,rep_fit,reps) = Seqs_freq_IAV(1);
    end

end
end

%% Export results -----------------------------------------------------------

mF1_IAV = mean(F1_IAV,3); % Mean mutant frequency after nG generations

% Save matrix as CSV
path = pwd;
date = sprintf(datestr(now,'mm-dd-yyyy_HH-MM-SS'));
dirName = ['Result_',date];
mkdir(path,dirName);
fullPath = fullfile(path,dirName,['Freq_Fix_',date,'.csv']);
writematrix(mF1_IAV,fullPath);

%% Plot heatmap -------------------------------------------------------------

figure
surf(W_factor, vMOI, mF1_IAV)
view(0,90)                     % View from top
colormap('jet')                % Color scheme
xlabel('Fitness')
ylabel('MOI')
colorbar
axis([W_factor(1) W_factor(end) vMOI(1) vMOI(end)])
ax = gca;
ax.FontSize = 20;
