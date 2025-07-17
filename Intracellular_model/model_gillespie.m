% This function performs one Gillespie step of the stochastic IAV replication model.
% It updates the system state vectors based on a randomly selected reaction, including special handling
% for virion production and protein interactions.

function [Si,SiM,SiMC,Frac_G,Frac_C,nv,Vir_time] = model_gillespie(Si,SiM,SiMC,Si1,SiM1,A1,dt,MR,MRM,IndS,Pos_syn,RNA_nuc,V0,nRNP,Rcrit,Frac_G,Frac_C,nv,Vir_time,Vir_Prod,NPs,PM_cyt,Prots,NumG,NumS,maxCycle,ti,Frac_I,MBind,M1,MParsM,Pol,NP,NEP,Reg,NumP,IndGen)

% Select next reaction based on Gillespie step
Ac = cumsum(A1(Rcrit))./sum(A1(Rcrit));
j1 = Rcrit(find(rand < Ac,1,'first'));

% Update global state vector Si with net stoichiometric change from the reaction
Si = Si + MR(:,j1)';

% Track number of reactions triggered (in the Gillespie step is always 1)
NumReacts = 1;

% Apply reaction to multi-variant state and cycle structure
[Si,SiM,SiMC] = ReactMOI(j1,A1,SiM1,IndS,Si,MRM,SiM,NumS,maxCycle,Pos_syn,SiMC,NumReacts,Vir_Prod,dt,RNA_nuc,V0,nRNP);

% If virion production occurred, update genotype tracking and virion info
if j1 == Vir_Prod
    [Si,SiM,SiMC,Frac_G,Frac_C,nv,Vir_time] = ReactVir(NumReacts,NPs,j1,A1,NumP,SiM,Si,PM_cyt,NumG,Frac_G,NumS,maxCycle,SiMC,Frac_C,Prots,IndGen,nv,dt,Vir_time,ti,Frac_I,Si1,SiM1);
end

% Check for negative concentrations (debugging)
if any(Si<0)
    Si
    Si1
    j1
    A1
end

% Re-synchronize Si and SiM in case of numerical drift
S = find(Si(1:NumS-1) ~= sum(SiM(:,1:NumS-1),1));
dS = Si(S) - sum(SiM(:,S),1);
Si(S) = Si(S) - dS;

% Apply protein binding reactions (M1 binding to RNPs)
for reac = 1:size(MBind,1)
    Rs = MBind(reac,:);
    [aP,bP] = find(MR(M1,Rs)<0);
    Prots_reac = (M1(unique(aP)))';

    if Si(Prots_reac) - sum(SiM(:,Prots_reac),1) ~= 0
        Params_reac = MParsM(:,Rs(1));
        SiM = ReactProts(Rs,A1,SiM1,SiM,dt,Prots_reac,Params_reac,Si);
    end
end

% Regulate polymerase, NP, and NEP proteins (e.g., binding to genomic complexes)
Prots_reac = Pol;
if Si(Prots_reac) - sum(SiM(:,Prots_reac),1) ~= 0
    Params_reac = Reg(3);
    SiM = ReactProtsGen(SiM1,SiM,Prots_reac,Params_reac,Si);
end

Prots_reac = NP;
if Si(Prots_reac) - sum(SiM(:,Prots_reac),1) ~= 0
    Params_reac = Reg(4);
    SiM = ReactProtsGen(SiM1,SiM,Prots_reac,Params_reac,Si);
end

Prots_reac = NEP;
if Si(Prots_reac) - sum(SiM(:,Prots_reac),1) ~= 0
    Params_reac = Reg(6);
    SiM = ReactProtsGen(SiM1,SiM,Prots_reac,Params_reac,Si);
end

end
