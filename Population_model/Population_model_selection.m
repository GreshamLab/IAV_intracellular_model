%% SELECTION SIMULATION IN IAV POPULATIONS
% -------------------------------------------------------------------------
% This script performs population-level simulations of Influenza A Virus (IAV)
% using pre-generated single-cell infection results (in this case, we're using
% the ones form Result_03-23-2025_08-14-13, but these can be generated with 
% the script Central_IAVGillespie_multiple). Unlike previous scripts,
% this one introduces a beneficial mutant variant with increased fitness. The
% simulation tracks its frequency over multiple viral generations under varying
% MOI and selection intensities.
% -------------------------------------------------------------------------

% =====================================================================
% This script is used to generate simulated data like the ones used 
% for the Figure 5A of the paper Segredo-Otero and Gresham 2025.
% =====================================================================

clear all
close all
clc

% Load simulation results
current = strcat(pwd,'/Result_03-23-2025_08-14-13');
path = current(end-19:end);
M = 20;
MOIs = 1:M;
NumG = 8;
Rep = 1000;
maxVir = 1e3;
Genomes_store = nan(NumG*Rep,M,M);
Fit_store = nan(Rep,M);

% Load pre-simulated single-cell data (genomes and fitness)
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

% Summarize frequencies and fitness across MOIs
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
    mFreq(1,m) = mean(Fit>10);
    vFreq(1,m) = sqrt(var(Fit>10));
    FreqAll(:,find(Fit == 0)) = [];
    Fit(find(Fit == 0)) = [];
    mFit(1,m) = mean(Fit);
    vFit(1,m) = sqrt(var(Fit));
    mFreqAll(:,m) = nanmean(FreqAll');
    vFreqAll(:,m) = sqrt(nanvar(FreqAll'));
end

%% SIMULATION PARAMETERS FOR MUTANT SELECTION
Total_differentes_seqs = 1e4;
Ini_title = 1e4;
Ini_freq = 1/Total_differentes_seqs;

Reps = 10;               % Number of simulation replicates per condition
Rep_fit = 2;             % Number of fitness levels tested
Rep_MOI = 2;             % Number of MOI levels tested
Fits = [1.5 5 1.5 5];    % Relative fitness of mutant in each scenario
MOIs = [1 1 20 20];      % Corresponding MOI values

Days = 14;               % Simulation duration (days)
nG = Days*2;             % Number of viral generations (2/day)
vMOI = linspace(1,20,Rep_fit);

F1_IAV = zeros(nG+1,Reps,Rep_fit*Rep_MOI);  % Stores mutant frequency over time
F1_IAV(1,:,:) = Ini_freq;

Segment = randi(8);      % Randomly choose one genome segment to simulate

% Iterate over fitness and MOI combinations
for rep_condition = 1:Rep_fit*Rep_MOI

    MOI = MOIs(rep_condition);
    Cells_per_generation = ones(1,nG).*round(1e5/MOI);
    w_factor = ones(1,Total_differentes_seqs); 
    w_factor(1) = Fits(rep_condition);          % Assign fitness advantage to mutant
    w_factor = w_factor./sum(w_factor);

    for reps = 1:Reps

        Seqs_freq = ones(1,Total_differentes_seqs).*Ini_freq;
        Num_seqs = mnrnd(Ini_title,Seqs_freq);
        Seqs_freq_IAV = Seqs_freq;

        for n = 1:nG

            Cells = Cells_per_generation(n);
            cMOI = MOI;
            Genomes = reshape(Genomes_store(:,1:cMOI,cMOI),[],cMOI,Rep);
            Infections_MOI = Genomes(Segment,:,:)./sum(Genomes(Segment,:,:),2);
            Infections_MOI(isnan(Infections_MOI)) = 0;

            Cells_population_IAV = randsample(Total_differentes_seqs,cMOI*Cells,true,Seqs_freq_IAV);
            Cells_population_IAV = reshape(Cells_population_IAV,cMOI,Cells);
            Rand_infs = randi(Rep,1,Cells);
            Infections_IAV = reshape(Infections_MOI(:,:,Rand_infs),cMOI,[]).*Fit_store(Rand_infs,cMOI)';

            Seqs_count_IAV = zeros(1,Total_differentes_seqs);
            for m = 1:cMOI
                Seqs_count_IAV(1,:) = Seqs_count_IAV(1,:) + ((accumarray(Cells_population_IAV(m,:)', Infections_IAV(m,:)', [Total_differentes_seqs, 1], @sum, 0))').*w_factor; 
            end

            Seqs_freq_IAV = Seqs_count_IAV./sum(Seqs_count_IAV); % Update frequencies
            F1_IAV(n+1,reps,rep_condition) = Seqs_freq_IAV(1);    % Store mutant frequency

        end
    end
end

%% PLOT MUTANT FREQUENCY OVER TIME -------------------------------------------
figure
hold on
plot(0:nG,F1_IAV(:,1,1)','--r','LineWidth',2)
plot(0:nG,F1_IAV(:,1,2)','-r','LineWidth',2)
plot(0:nG,F1_IAV(:,1,3)','--b','LineWidth',2)
plot(0:nG,F1_IAV(:,1,4)','-b','LineWidth',2)
plot(0:nG,F1_IAV(:,:,1)','--r','LineWidth',2)
plot(0:nG,F1_IAV(:,:,2)','-r','LineWidth',2)
plot(0:nG,F1_IAV(:,:,3)','--b','LineWidth',2)
plot(0:nG,F1_IAV(:,:,4)','-b','LineWidth',2)
set(gca,'YScale','log','Ylim',[1e-6 1],'Xlim',[0 nG])
xlabel('Time (generations)')
ylabel('HFV frequency')
ax = gca;
ax.FontSize = 15;
legend('w_i = 1.5, MOI = 1','w_i = 5, MOI = 1', ...
       'w_i = 1.5, MOI = 20','w_i = 5, MOI = 20')
