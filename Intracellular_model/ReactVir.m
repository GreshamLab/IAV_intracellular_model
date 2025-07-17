% This function models the production of a single virion within the stochastic IAV replication model.
% It samples proteins and genomic segments from cytoplasmic pools (Si and SiM), and stores
% the virion genotype and replication cycle composition (Frac_G, Frac_C), tracking information (nv, Vir_time), and
% cycle-resolved matrix (SiMC). Ensures stoichiometric consistency and prevents negative concentrations.

function [Si,SiM,SiMC,Frac_G,Frac_C,nv,Vir_time] = ReactVir(NumReacts,NPs,j1,A1,NumP,SiM,Si,PM_cyt,NumG,Frac_G,NumS,maxCycle,SiMC,Frac_C,Prots,IndGen,nv,dt,Vir_time,ti,Frac_I,Si1,SiM1)

% Total number of virions to be produced in this call
newV = NumReacts;

% Compute total protein copies required across all new virions
TotalProt = NPs(2,:) .* newV;

% Estimate how proteins are distributed across MOIs
fracProt = 1 - exp(-A1(j1) * dt * SiM(:,Prots) ./ sum(SiM(:,Prots), 1));
if sum(fracProt) > 0.1
    fracProt = SiM(:,Prots) ./ sum(SiM(:,Prots), 1);
end
fracProt = fracProt ./ sum(fracProt, 1);

% Sample how many proteins are drawn from each MOI
rProts = (mnrnd(TotalProt', fracProt', NumP))';
rProts(isnan(rProts)) = 0;

% Subtract sampled proteins from each MOI (SiM) and global pool (Si)
if length(SiM(:,Prots)) == length(rProts)
    SiM(:,Prots) = SiM(:,Prots) - rProts;
end
Si(Prots) = Si(Prots) - TotalProt;

% Enforce non-negativity: correct any negative values due to stochastic sampling
count = 0;
while any(SiM(:) < 0)
    count = count + 1;
    if count > 100
        break
    end

    NegP = find(Si(Prots) == 0);
    SiM(:,Prots(NegP)) = 0;

    [Mneg, Sneg] = find(SiM < 0);
    for sN = 1:length(Sneg)
        negs = find(SiM(:,Sneg(sN)) < 0);
        SiM(negs, Sneg(sN)) = 0;
        realSM = Si(Sneg(sN));
        if sum(SiM(:,Sneg(sN))) > 0
            SiM(:,Sneg(sN)) = realSM * (SiM1(:,Sneg(sN)) ./ sum(SiM1(:,Sneg(sN))));
        else
            SiM(:,Sneg(sN)) = 0;
        end
    end
end

% Re-add proteins to global pool (Si) to maintain mass balance
Si(Prots) = Si(Prots) + TotalProt;

% Select which genome segments to include in each new virion
Av = sum(Si1(PM_cyt) > 0);
Num = fliplr(1:8);

for n = 1:newV
    % Determine number of genome segments for this virion
    Number_segments = mnrnd(1, Frac_I);
    Segs = Num(find(Number_segments));
    Segs(Segs > Av) = Av;

    % Sample which genome segments will be packaged from available cytoplasmic pools
    sPM = randsample_WithoutReplacement(1:NumG, Segs, Si1(PM_cyt) ./ sum(Si1(PM_cyt)));
    ssPM = size(sPM);
    Si(PM_cyt(sPM)) = Si(PM_cyt(sPM)) - 1;

    % Select donor variant (MOI) for each segment based on availability
    Genomes_vir = floor(SiM(:,PM_cyt(sPM)));
    Genomes_vir = Genomes_vir ./ sum(Genomes_vir, 1);
    Genomes_vir = cumsum(Genomes_vir, 1);
    j_gens = sum(rand(1,ssPM(2)) > Genomes_vir) + 1;

    % Remove selected genome segments from each MOI
    Genomes_vir = SiM(:,PM_cyt(sPM));
    Genomes_vir(IndGen(:,sPM) == j_gens) = Genomes_vir(IndGen(:,sPM) == j_gens) - 1;
    SiM(:,PM_cyt(sPM)) = Genomes_vir;
    if sum(SiM(:) < 0) > 0
        break
    end

    % Record genotype composition of the virion
    nv = nv + 1;
    Frac_G(sPM, nv) = j_gens;

    % Determine cycle of replication for each genome segment
    for j = 1:ssPM(2)
        gens = sPM(j);
        SiMC2 = reshape(SiMC(j_gens(j),:), NumS, maxCycle)';
        Genomes_cycle = SiMC2(:,PM_cyt(gens));

        % Redistribute if negative values appear in the cycle-specific counts
        count = 0;
        while any(Genomes_cycle < 0)
            count = count + 1;
            if count > 100
                break
            end
            c1 = find(Genomes_cycle < 0);
            others = setdiff(1:maxCycle, c1);
            P_M = Genomes_cycle(others) ./ sum(Genomes_cycle(others));
            P_M = P_M ./ sum(P_M);
            if ~any(isnan(P_M))
                Genomes_cycle(others) = Genomes_cycle(others) + sum(Genomes_cycle(c1)) * P_M;
                Genomes_cycle(c1) = 0;
            end
        end

        % Sample cycle for the segment being incorporated
        P_M = 1 - exp(-A1(j1) * dt * Genomes_cycle ./ sum(Genomes_cycle));
        P_M = P_M ./ sum(P_M);
        j3 = mnrnd(1, P_M);
        if ~any(isnan(j3))
            Frac_C(gens, nv) = find(j3);
            Genomes_cycle = Genomes_cycle - j3';
        end

        % Final redistribution if negative values still exist
        while any(Genomes_cycle < 0)
            c1 = find(Genomes_cycle < 0);
            others = setdiff(1:maxCycle, c1);
            P_M = Genomes_cycle(others) ./ sum(Genomes_cycle(others));
            P_M = P_M ./ sum(P_M);
            if ~any(isnan(P_M))
                Genomes_cycle(others) = Genomes_cycle(others) + sum(Genomes_cycle(c1)) * P_M;
                Genomes_cycle(c1) = 0;
            end
        end

        % Update SiMC with adjusted values
        SiMC2(:,PM_cyt(gens)) = Genomes_cycle;
        if any(any(SiMC2(:, PM_cyt) < 0))  || any(any(isnan(SiMC2(:, PM_cyt))))
            break
        end
        SiMC(j_gens(j), :) = reshape(SiMC2', 1, NumS * maxCycle);
    end

    % Track virion production time
    Vir_time(1,nv) = ti;

    % Debug: verify integrity of SiMC vs SiM
    MOI = size(SiM, 1);
    for m = 1:MOI
        SiMC2 = reshape(SiMC(m,:), NumS, maxCycle)';
        if prod(round(sum(SiMC2(:,1:35))) == round(SiM(m,1:35))) == 0
            find((sum(SiMC2) == SiM(m,:)) == 0)
            m
            break
        end
    end
end

end
