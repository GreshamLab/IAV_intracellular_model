% This function updates the replication cycle-resolved state matrix (SiMC) for each MOI,
% distributing the effect of a reaction (j1) across replication cycles based on probabilities 
% derived from their current segment content. It ensures stoichiometric consistency and avoids 
% negative concentrations. If the reaction is the synthesis of pRNPs, the cycle of the produce is +1.

function SiMC = RepCycles(j2, NumS, maxCycle, A1, j1, IndS, Pos_syn, SiMC, cS, MRM, dt, SiM)

% Get the indices of non-zero reactions to process
IndJ = find(j2);
j2 = j2(IndJ);

% Loop over each reacting MOI
for sj = 1:size(j2',1)
    
    % Extract replication cycle matrix for current MOI (SiMC2 is maxCycle x NumS)
    SiMC2 = reshape(SiMC(IndJ(sj),:), NumS, maxCycle)';

    % Compute probability distribution over cycles
    P_M = 1 - exp(-A1(j1) * dt * SiMC2(:,IndS(1,j1)) ./ sum(SiMC2(:,IndS(1,j1))));
    if sum(P_M) > 0.1
        P_M = SiMC2(:,IndS(1,j1)) ./ sum(SiMC2(:,IndS(1,j1)));
    end
    P_M = P_M ./ sum(P_M);

    % Separate integer and decimal part of j2
    j2I = floor(j2(sj));
    j2d = j2(sj) - j2I;

    % Multinomial allocation of reaction over cycles
    j3 = mnrnd(j2I, P_M)';

    % Increasing the cycle number of the reaction is a synthesis of pRNPs
    if ismember(j1, Pos_syn)
        j3(2:end) = j3(1:end-1);
        j3(1) = 0;
    end

    % Apply update to SiMC2 (cycle-resolved matrix)
    if ~any(isnan(j3))
        JS = size((MRM(cS, j1) .* j3')');
        if JS(2) ~= maxCycle
            SiMC2(:,cS) = SiMC2(:,cS) + (MRM(cS, j1) * j3')';
        end
    end

    % Check and correct for negative values
    while sum(any(SiMC2(:,cS) < 0)) > 0
        [aN, bN] = find(SiMC2(:,cS) < 0);
        diffC = sum(abs(SiMC2(aN, cS(1))));
        if length(cS) > 1
            SiMC2(aN, cS(2)) = SiMC2(aN, cS(2)) + SiMC2(aN, cS(1));
        end
        SiMC2(aN, cS(1)) = 0;

        % Re-distribute reaction proportionally
        if sum(SiMC2(:,IndS(1,j1))) > 0
            P_M = 1 - exp(-A1(j1) * dt * SiMC2(:,IndS(1,j1)) ./ sum(SiMC2(:,IndS(1,j1))));
            if sum(P_M) > 0.1
                P_M = SiMC2(:,IndS(1,j1)) ./ sum(SiMC2(:,IndS(1,j1)));
            end
            P_M = P_M ./ sum(P_M);
            if ~any(isnan(P_M))
                SiMC2(:,cS) = SiMC2(:,cS) + diffC .* (MRM(cS,j1) * P_M')';
            end
        end
    end

    if sum(SiMC(:,cS) < 0) > 0
        break
    end

    % Handle remaining decimal part of j2
    if j2d > 0
        j3 = mnrnd(1, P_M)';

        if ismember(j1, Pos_syn)
            j3(2:end) = j3(1:end-1);
            j3(1) = 0;
        end

        if ~any(isnan(j3))
            JS = size((MRM(cS, j1) .* j3')');
            if JS(2) ~= maxCycle
                SiMC2(:,cS) = SiMC2(:,cS) + j2d .* (MRM(cS, j1) * j3')';
            end
        end
    end

    % Repeat correction for potential negative values after decimal update
    while sum(any(SiMC2(:,cS) < 0)) > 0
        [aN, bN] = find(SiMC2(:,cS) < 0);
        diffC = sum(abs(SiMC2(aN, cS(1))));
        if length(cS) > 1
            SiMC2(aN, cS(2)) = SiMC2(aN, cS(2)) + SiMC2(aN, cS(1));
        end
        SiMC2(aN, cS(1)) = 0;

        if sum(SiMC2(:,IndS(1,j1))) > 0
            P_M = 1 - exp(-A1(j1) * dt * SiMC2(:,IndS(1,j1)) ./ sum(SiMC2(:,IndS(1,j1))));
            if sum(P_M) > 0.1
                P_M = SiMC2(:,IndS(1,j1)) ./ sum(SiMC2(:,IndS(1,j1)));
            end
            P_M = P_M ./ sum(P_M);
            if ~any(isnan(P_M))
                SiMC2(:,cS) = SiMC2(:,cS) + diffC .* (MRM(cS,j1) * P_M')';
            end
        end
    end

    if sum(SiMC(:) < 0) > 0
        break
    end

    if sum(isnan(SiMC2(:)))
        break
    end

    % Save updated cycle-specific state back into SiMC
    SiMC(IndJ(sj), :) = reshape(SiMC2', 1, NumS * maxCycle);

end

end
