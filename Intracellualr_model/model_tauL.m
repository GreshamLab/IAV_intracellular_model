% This function performs one tau-leaping step in the stochastic IAV replication model.
% It updates system states based on Poisson-distributed firing of non-critical reactions,
% including virion production and protein interactions.

function [Si,SiM,SiMC,Frac_G,Frac_C,nv,Vir_time] = model_tauL(Si,SiM,SiMC,Si1,SiM1,A1,dt,MR,MRM,IndS,Pos_syn,RNA_nuc,V0,nRNP,Rnocrit,Frac_G,Frac_C,nv,Vir_time,Vir_Prod,NPs,PM_cyt,Prots,NumG,NumS,maxCycle,ti,Frac_I,MBind,M1,MParsM,Pol,NP,NEP,Reg,NumP,IndGen)

% Sample the number of firings for each non-critical reaction using Poisson distribution
Nreacts = poissrnd(A1(Rnocrit) .* dt);
Nr = find(Nreacts);  % Indices of reactions that fired

% If no reactions fire, return immediately
if sum(Nreacts) == 0
    return
end

% Iterate over each reaction that fired
for sr = 1:length(Nr)
    j1 = Rnocrit(Nr(sr));           % Get reaction index
    NumReacts = Nreacts(Nr(sr));    % Number of firings for this reaction

    % Apply the reaction to multi-variant (MOI) and cycle-specific state
    [Si,SiM,SiMC] = ReactMOI(j1,A1,SiM1,IndS,Si,MRM,SiM,NumS,maxCycle,Pos_syn,SiMC,NumReacts,Vir_Prod,dt,RNA_nuc,V0,nRNP);

    % If virion production occurred, trigger genotype tracking and virion updates
    if j1 == Vir_Prod
        VirP = find(Rnocrit == Vir_Prod);
        maxV = floor(min(max(SiM(:,PM_cyt))));  % Maximum virions constrained by available RNPs

        % Limit virion production if insufficient material
        if Nreacts(VirP) >= maxV
            Nreacts(VirP) = maxV - 1;
        end

        NumReacts = Nreacts(Nr(sr));
        [Si,SiM,SiMC,Frac_G,Frac_C,nv,Vir_time] = ...
            ReactVir(NumReacts,NPs,j1,A1,NumP,SiM,Si,PM_cyt,NumG, ...
                     Frac_G,NumS,maxCycle,SiMC,Frac_C,Prots, ...
                     IndGen,nv,dt,Vir_time,ti,Frac_I,Si1,SiM1);
    end
end

% Apply the global state update from all Poisson-sampled reactions
Si = Si + sum(MR(:,Rnocrit(Nr)) .* Nreacts(Nr), 2)';

% Re-synchronize Si with SiM in case of numeric drift or inconsistency
S = find(Si(1:NumS-1) ~= sum(SiM(:,1:NumS-1),1));
dS = Si(S) - sum(SiM(:,S),1);
Si(S) = Si(S) - dS;

% Debugging: check for negative states
if any(Si < 0)
    Nreacts
    Rnocrit(Nr)
    Si
    Si1
end

% Apply protein binding reactions (M1 binding to RNPs)
for reac = 1:size(MBind,1)
    Rs = MBind(reac,:);
    [aP,~] = find(MR(M1,Rs) < 0);
    Prots_reac = (M1(unique(aP)))';

    if Si(Prots_reac) - sum(SiM(:,Prots_reac),1) < 0
        Params_reac = MParsM(:,Rs(1));
        SiM = ReactProts(Rs,A1,SiM1,SiM,dt,Prots_reac,Params_reac,Si);
    end
end

% Regulate polymerase, NP, and NEP proteins if needed
Prots_reac = Pol;
if Si(Prots_reac) - sum(SiM(:,Prots_reac),1) < 0
    Params_reac = Reg(3);
    SiM = ReactProtsGen(SiM1,SiM,Prots_reac,Params_reac,Si);
end

Prots_reac = NP;
if Si(Prots_reac) - sum(SiM(:,Prots_reac),1) < 0
    Params_reac = Reg(4);
    SiM = ReactProtsGen(SiM1,SiM,Prots_reac,Params_reac,Si);
end

Prots_reac = NEP;
if Si(Prots_reac) - sum(SiM(:,Prots_reac),1) < 0
    Params_reac = Reg(6);
    SiM = ReactProtsGen(SiM1,SiM,Prots_reac,Params_reac,Si);
end

end
