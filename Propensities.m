function A1 = Propensities(Si1, IndS, MPars, Vir_Prod, Params, PM_cyt, Prots, NPs, ...
                           Translations, mRNA, Transcriptions, Pol, ...
                           Neg_syn, Pos_syn, Reg, NP, MR, Si, MBind, NEP, SiM1, NumS)
% PROPENSITIES computes the reaction propensities A1 for all reactions at the current time step.
% Includes mass-action terms and analytical approximations for fast-switching intermediates
% (e.g., naked RNA, RNA-Pol complexes, and RNP-M1 intermediates) based on QSSA assumptions.

    % -------------------------------
    % General mass-action propensities
    % -------------------------------
    A1 = prod([Si1(1,IndS(1,:)); Si1(1,IndS(2,:))], 1) .* MPars;

    % -------------------------------
    % Virion production
    % -------------------------------
    A1(Vir_Prod) = Params(11) * ...
                   min(Si1(PM_cyt(:))) * ...
                   prod(Si1(Prots) ./ (Si1(Prots) + Params(23) * NPs(1,:)));

    % Disable virion production if no genotype has a full set of PM_cyt segments
    if prod(sum(SiM1(:, PM_cyt) >= 1) > 0) == 0
        A1(Vir_Prod) = 0;
    end

    % -------------------------------
    % Eliminate reactions using extinct species
    % -------------------------------
    B = find(sum(SiM1(:, 1:NumS-1) >= 1, 1) == 0);
    for b = 1:length(B)
        [~, R] = find(B(b) == IndS);
        A1(R) = 0;
    end

    % -------------------------------
    % Translation (special case for segment 1)
    % -------------------------------
    A1(Translations(1)) = min(Si1(mRNA(1, 1:3)))' * MPars(Translations(1));

    % -------------------------------
    % Transcription feedback approximation
    % -------------------------------
    % Simple negative feedback from polymerase accumulation
    A1(Transcriptions) = A1(Transcriptions) ./ (1 + Si1(Pol) / Reg(7));

    % -------------------------------
    % QSSA-based approximation for RNP formation (from naked RNA)
    % -------------------------------
    % These reactions simulate the net formation of RNPs (both positive and negative strand)
    % via rapid intermediate states: naked RNA → RNA-Pol → RNP
    % Each step has an associated survival probability under high degradation

    A1(Neg_syn) = A1(Neg_syn) ./ ...
        (1 + Reg(1) / (Reg(3) * (Si1(Pol) + 1e-10))) ./ ...
        (1 + Reg(2) / (Reg(4) * (Si1(NP) + 1e-10)));

    A1(Pos_syn) = A1(Pos_syn) ./ ...
        (1 + Reg(1) / (Reg(3) * (Si1(Pol) + 1e-10))) ./ ...
        (1 + Reg(2) / (Reg(4) * (Si1(NP) + 1e-10)));

    % These approximations represent:
    % PRNA  = k_bind_pol / (k_bind_pol + k_deg_nakedRNA)
    % PRNAP = k_bind_NP / (k_bind_NP + k_deg_RNA-pol)
    % The overall RNP formation propensity is the original one multiplied by PRNA * PRNAP

    % -------------------------------
    % QSSA-based approximation for nuclear export of RNPs
    % -------------------------------
    % RNP-M1 intermediates can be degraded or exported after binding NEP.
    % This term approximates the probability of successful export:
    % PRNPM1 = k_bind_NEP / (k_bind_NEP + k_deg_RNP-M1)

    A1(MBind) = A1(MBind) ./ (1 + Reg(5) / (Reg(6) * (Si1(NEP) + 1e-10)));

    % -------------------------------
    % Physical feasibility check
    % -------------------------------
    Mcrit = Si' + MR;
    [~, bC] = find(Mcrit < 0);
    Rnull = unique(bC);
    A1(Rnull) = 0;

end

                                                                                                                 