function [dt, dt_noC, dt_C, Rcrit, Rnocrit] = dtCal(Si, PM_cyt, MRC, Vir_Prod, NumR, MS, A1A, MR, e, gi)
% Calculates time steps for stochastic simulation using a hybrid method:
% - dt: global Gillespie step
% - dt_noC: tau-leap step for non-critical reactions
% - dt_C: Gillespie step for critical reactions
% - Rcrit: critical reactions (risk species depletion)
% - Rnocrit: non-critical reactions (safe for tau-leaping)

% Identify reactions that would deplete species
sPM = find(Si(PM_cyt));
MRC(PM_cyt(sPM), Vir_Prod) = -1;
MRC(MRC > 0) = 0;

Mcrit = Si' + MRC * 11;
[~, bC] = find(Mcrit < 0);
Rcrit = unique(bC);
sR1 = 1:NumR;
Rnocrit = setdiff(sR1, Rcrit);

% Compute expected drift (mu) and variance (sigma) for non-critical reactions
dM1 = sum((A1A(Rnocrit) .* MR(:, Rnocrit))', 1);
mu = max(e .* Si ./ gi, 1) ./ abs(dM1);

dM2 = sum((A1A(Rnocrit) .* MR(:, Rnocrit).^2)', 1);
sig = max(e .* Si ./ gi, 1).^2 ./ dM2;

% Time steps
dt = exprnd(1 / sum(A1A));             % Gillespie for all
dt_noC = min([mu, sig]);              % Safe tau step for non-critical
if dt_noC > 10*60; dt_noC = 0; end     % Cap large tau
dt_C = exprnd(1 / sum(A1A(Rcrit)));    % Gillespie step for critical reactions

end
