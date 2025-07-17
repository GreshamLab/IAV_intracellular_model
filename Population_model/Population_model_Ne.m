%% SINGLE-GENERATION VIRAL INFECTION SIMULATION
% -------------------------------------------------------------------------
% This script simulates a single round of infection to estimate the effective
% population size (Ne) of viral populations under two models:
% 1. Null model: all virions are equally fit and all infections are successful.
% 2. IAV model: uses simulated genome and fitness data from prior 
% single-cell runs (in this case, we're using the ones form 
% Result_03-23-2025_08-14-13, but these can be generated with the script 
% Central_IAVGillespie_multiple). Each simulation tracks 10^4 neutral 
% variants (barcodes) and computes Ne by comparing frequencies before and after infection.
% -------------------------------------------------------------------------


% ==============================================================================
% This script is used to generate simulated data for the Figure 4A of the paper 
% Segredo-Otero and Gresham 2025.
% ==============================================================================

clear all
close all
clc

%% LOAD SIMULATION DATA ---------------------------------------------------
current = strcat(pwd,'/Result_03-23-2025_08-14-13');
path = current(end-19:end);
M = 20;                          % Number of different MOI conditions
MOIs = 1:M;
NumG = 8;                        % Number of genome segments
Rep = 1000;                      % Number of replicates
Genomes_store = nan(NumG*Rep,M,M);
Fit_store = nan(Rep,M);

% Load pre-generated data for each MOI
for m = 1:M
    Genomes_name = strcat(current,'/MOI_',num2str(MOIs(m)),'_Genomes',path,'.csv');
    Fit_name = strcat(current,'/MOI_',num2str(MOIs(m)),'_Fit',path,'.csv');
    
    Genomes = csvread(Genomes_name);
    Genomes_store(:,1:MOIs(m),MOIs(m)) = Genomes;

    Fit = csvread(Fit_name);
    Fit = Fit(1,:);
    Fit(isnan(Fit)) = 0; Fit(Fit == Inf) = 0;
    Fit_store(:,MOIs(m)) = Fit;
end

%% COMPUTE AVERAGE FITNESS AND FREQUENCY STATS ----------------------------
mFreqAll = zeros(M,M);
mFit = zeros(1,M);
vFreqAll = zeros(M,M);
vFit = zeros(1,M);
mFreq = zeros(1,M);
vFreq = zeros(1,M);

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

    mFreq(1,m) = mean(Fit>50);
    vFreq(1,m) = sqrt(var(Fit>50));
    FreqAll(:,Fit == 0) = [];
    Fit(Fit == 0) = [];

    mFit(1,m) = mean(Fit);
    vFit(1,m) = sqrt(var(Fit));

    mFreqAll(:,m) = nanmean(FreqAll');
    vFreqAll(:,m) = sqrt(nanvar(FreqAll'));
end

%% SIMULATE SINGLE ROUND OF INFECTION -------------------------------------
Total_differentes_seqs = 1e4;
Ini_freq = 1/Total_differentes_seqs;
Seqs_freq = ones(1,Total_differentes_seqs).*Ini_freq;  % Equal initial freq

Rep_MOI = M;
Reps = 1000;  % Number of simulated replicates
MOI = linspace(1,M,Rep_MOI);
Cells_per_generation = 1e4;

NT_vir = zeros(Reps,Rep_MOI);
NT_IAV = zeros(Reps,Rep_MOI);
NT_real = zeros(Reps,Rep_MOI);

for rep_MOI = 1:Rep_MOI
    cMOI = MOI(rep_MOI);
    Cells = Cells_per_generation;
    Genomes = reshape(Genomes_store(:,1:cMOI,cMOI),[],cMOI,Rep);
    Infections_MOI = Genomes(6,:,:)./sum(Genomes(6,:,:),2);
    Infections_MOI(isnan(Infections_MOI)) = 0;

    for rep = 1:Reps
        Cells_population = randsample(Total_differentes_seqs,cMOI*Cells,true,Seqs_freq);
        Cells_population = reshape(Cells_population,cMOI,Cells);
        Cells_population_vir = Cells_population;
        Cells_population_IAV = Cells_population;
        Rand_infs = randi(Rep,1,Cells);
        Infections_IAV = reshape(Infections_MOI(:,:,Rand_infs),cMOI,[]).*Fit_store(Rand_infs,cMOI)';

        Seqs_count_vir = zeros(1,Total_differentes_seqs);
        Seqs_count_IAV = zeros(1,Total_differentes_seqs);

        for m = 1:cMOI
            Seqs_count_IAV(1,:) = Seqs_count_IAV(1,:) + (accumarray(Cells_population_IAV(m,:)', Infections_IAV(m,:)', [Total_differentes_seqs, 1], @sum, 0))';
            Seqs_count_vir(1,:) = Seqs_count_vir(1,:) + histcounts(Cells_population_vir(m,:), 1:Total_differentes_seqs + 1);
        end

        Seqs_freq_vir = Seqs_count_vir./sum(Seqs_count_vir);
        Seqs_freq_IAV = Seqs_count_IAV./sum(Seqs_count_IAV);
        if Seqs_count_IAV == 0; Seqs_freq_IAV(:) = 0; end 
        if Seqs_count_vir == 0; Seqs_freq_vir(:) = 0; end

        NT_vir(rep,rep_MOI) = mean((Seqs_freq.*(1-Seqs_freq)))./mean((Seqs_freq_vir - Seqs_freq).^2);
        NT_IAV(rep,rep_MOI) = mean((Seqs_freq.*(1-Seqs_freq)))./mean((Seqs_freq_IAV - Seqs_freq).^2);
        NT_real(rep,rep_MOI) = Cells*cMOI;
    end
end

%% CALCULATE STATISTICS ---------------------------------------------------
mNT_vir = mean(NT_vir,1);
vNT_vir = sqrt(var(NT_vir,0,1));
mNT_IAV = mean(NT_IAV,1);
vNT_IAV = sqrt(var(NT_IAV,0,1));
mNT_real = mean(NT_real,1);

Frac_NT = NT_IAV./NT_real;
mFrac = mean(Frac_NT,1);
vFrac = sqrt(var(Frac_NT,0,1));

%% PLOT: Ne COMPARISON FOR IAV VS NULL MODEL ------------------------------
figure
errorbar(MOI,mNT_vir,vNT_vir,'-*r','DisplayName','Null model')
hold on
errorbar(MOI,mNT_IAV,vNT_IAV,'-*b','DisplayName','IAV model')
plot(MOI,mNT_vir,'-*r',MOI,mNT_IAV,'-*b','LineWidth',2,'HandleVisibility','off')
hold off
set(gca,'YScale','log')
xlabel('MOI')
ylabel('\it{N_e}')
ax = gca;
ax.FontSize = 20;
legend('show', 'Location', 'best')
