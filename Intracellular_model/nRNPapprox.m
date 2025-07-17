function [Si,SiM,SiMC,A1,dt,dt_noC,dt_C,Rcrit,Rnocrit] = nRNPapprox(Si,SiM,SiMC,A1,MR,SiM1,nRNP,PM_cyt,Params,NumS,maxCycle,IndS,MBind,Neg_syn,Prots,M1,Pol,NP,NEP,Reg,Vir_Prod,NumR,MS,e,gi,MOI,MParsM,NPs,dt,dt_noC,dt_C)
% NRNPPROX  Quasi‑steady‑state (QSSA) shortcut for nRNP intermediates.
%
% This routine is called when nRNP‑related reactions (replication + M1 binding)
% limit the global step size (dt).  Instead of simulating naked RNA → nRNP and
% nRNP → nRNP‑M1 explicitly, we:
%   1) Estimate a larger safe timestep (dt_new) without those reactions.
%   2) Compute nRNP abundance deterministically as
%        nRNP = k_rep / (k_bind+ k_deg).
%   3) Propagate stoichiometric changes to cytoplasmic pools & proteins.
%   4) Redistribute updated counts across MOI variants and replication cycles
%      in proportion to their previous shares.
%   5) Restore full propensities and updated dt values.
%
% Inputs and outputs match the main simulation state vectors.
% -------------------------------------------------------------------------

%% 1.  Remove nRNP‑related reactions from propensity estimation
A1A = A1;                          % Save full propensities
A1A([Neg_syn MBind]) = 0;          % Mask replication & M1‑binding
[~,ReactsnRNP] = find(MR(nRNP,:) < 0); % Any other nRNP consumers
A1A(ReactsnRNP) = 0;

%% 2.  Recalculate tentative time steps without those reactions
dt_old = max([dt_noC dt*10]);
[dt1,dt_noC1,dt_C1,Rcrit,Rnocrit] = dtCal(Si,PM_cyt,MR,Vir_Prod, ...
    NumR,MS,A1A,MR,e,gi);

% Choose conservative new step
if dt_noC1 < dt1*10
    dt_new = dt1;                 % Gillespie dominates
else
    dt_new = min([dt_noC1 dt_C1]);% Hybrid tau‑leap dominates
end

%% 3.  Apply QSSA only if we gain >3× step size
if dt_new > 3*dt_old

    % Update global dt values used by caller
    dt      = dt1;
    dt_noC  = dt_noC1;
    dt_C    = dt_C1;

    %% 3a.  Deterministic nRNP concentration
    Si(nRNP) = A1(Neg_syn) ./ (A1(MBind)./Si(nRNP) + Params(27));
    Si(isnan(Si)) = 0;            % Clean numerical noise

    %% 3b.  Net stoichiometric change from masked reactions over dt_new
    dM1 = (A1([Neg_syn MBind]).*MR(:,[Neg_syn MBind]))';
    Si(PM_cyt) = Si(PM_cyt) + sum(dM1(:,PM_cyt),1).*dt_new;
    Si(Prots)  = Si(Prots)  + sum(dM1(:,Prots),1).*dt_new;
    Si(Si < 0) = 0;

    %% 3c.  Redistribute nRNP counts across MOI variants
    SiM(:,nRNP) = Si(nRNP).*SiM(:,nRNP)./sum(SiM(:,nRNP),1);
    SiM(:,nRNP(A1(Neg_syn)==0)) = 0;          % Zero if no replication
    SiM(:,nRNP(Si(nRNP)==0))    = 0;          % Zero if none left

    %% 3d.  Update PM_cyt per variant using weighted probabilities
    P_M = 1 - exp(-A1.*dt_new.*SiM1(:,IndS(1,:))./sum(SiM1(:,IndS(1,:))));
    P_M = sum(P_M,1);                  % Helper, avoids NaNs
    dSiM = P_M(:,[Neg_syn MBind])*dM1.*dt_new;
    SiM(:,PM_cyt) = SiM(:,PM_cyt) + dSiM(:,PM_cyt);

    %% 3e.  Per‑cycle redistribution in SiMC
    for m = 1:MOI
        SiMC2         = ones(maxCycle,NumS+1);
        SiMC2(:,1:NumS) = reshape(SiMC(m,:),NumS,maxCycle)';
        P_C = 1 - exp(-A1.*dt.*SiMC2(:,IndS(1,:))./sum(SiMC2(:,IndS(1,:))));
        P_C = sum(P_C,1);
        dM1c = (A1([Neg_syn MBind]).*MR(:,[Neg_syn MBind]))'.*P_M(m,[Neg_syn MBind])';
        dSiMC = P_C(:,[Neg_syn MBind])*dM1c*dt_new;
        SiMC2(:,PM_cyt) = SiMC2(:,PM_cyt) + dSiMC(:,PM_cyt);

        % Keep genotype proportions for nRNP across cycles
        nRNPp = find(SiM(m,nRNP));
        if ~isempty(nRNPp)
            SiMC2(:,nRNP(nRNPp)) = SiM(m,nRNP(nRNPp)).*SiMC2(:,nRNP(nRNPp))./ ...
                                   sum(SiMC2(:,nRNP(nRNPp)),1);
        end
        SiMC(m,:) = reshape(SiMC2(:,1:NumS)',1,NumS*maxCycle);
    end

    %% 3f.  Protein balancing.

    [aP,bP] = find(MR(M1,MBind)<0);
    Prots_reac = (M1(unique(aP)))';
    if Si(Prots_reac) - sum(SiM(:,Prots_reac),1) < 0
        Params_reac = MParsM(:,MBind(1));
        SiM = ReactProtsDet(MBind,A1,SiM,dt_new,Prots_reac,Params_reac,Si);
    end
    
    Prots_reac = Pol;
    Params_reac = Reg(3);
    if Si(Prots_reac) - sum(SiM(:,Prots_reac),1) < 0
        SiM = ReactProtsGenDet(SiM,Prots_reac,Params_reac,Si);
    end
    
    Prots_reac = NP;
    Params_reac = Reg(4);
    if Si(Prots_reac) - sum(SiM(:,Prots_reac),1) < 0
        SiM = ReactProtsGenDet(SiM,Prots_reac,Params_reac,Si);
    end
    
    Prots_reac = NEP;
    Params_reac = Reg(6);
    if Si(Prots_reac) - sum(SiM(:,Prots_reac),1) < 0
        SiM = ReactProtsGenDet(SiM,Prots_reac,Params_reac,Si);
    end
    
    S = find(Si(1:NumS-1) ~= sum(SiM(:,1:NumS-1),1));
    dS = Si(S) - sum(SiM(:,S),1);
    Si(S) = Si(S) - dS;
    
    %% 4.  Restore full propensities and update virion‑production term
    A1 = A1A;                           % Unmask
    Si1 = [Si 1];                       % Pad for indexing safety
    A1(Vir_Prod) = Params(11).*min(Si(PM_cyt)).* ...
                   prod(Si(Prots)./(Si(Prots)+Params(23).*NPs(1,:)));
    if ~any(all(SiM1(:,PM_cyt) >= 1,2))
        A1(Vir_Prod) = 0;               % No complete genome sets
    end
    
    %% 5.  Re‑evaluate dt ranges with updated state
    [dt1,dt_noC1,dt_C1,Rcrit,Rnocrit] = dtCal(Si,PM_cyt,MR,Vir_Prod, ...
        NumR,MS,A1,MR,e,gi);
    if sum(A1(Rcrit))   == 0, dt_C   = Inf; end
    if sum(A1(Rnocrit)) == 0, dt_noC = Inf; end

end

end
