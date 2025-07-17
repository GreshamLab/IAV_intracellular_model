% This function applies the effect of a given reaction (j1) on the MOI-specific system state.
% It handles reaction distribution across MOis (SiM), and
% updates replication cycles (SiMC) and nuclear RNA import. This is used within both Gillespie
% and tau-leaping steps of the IAV model.

function [Si,SiM,SiMC] = ReactMOI(j1,A1,SiM1,IndS,Si,MRM,SiM,NumS,maxCycle,Pos_syn,SiMC,NumReacts,Vir_Prod,dt,RNA_nuc,V0,nRNP)

% Skip if the reaction is virion production (handled elsewhere)
if j1 ~= Vir_Prod
    % Estimate relative contribution of each MOI index to the reactants involved
    ReacsS = floor(SiM1(:,IndS(1,j1))) ./ sum(floor(SiM1(:,IndS(1,j1))));
    
    % Compute probabilities of reaction firing per MOI using propensities
    P_M = 1 - exp(-A1(j1) * dt * ReacsS);
    if sum(P_M) > 0.1
        P_M = ReacsS;  % fallback to distribution if exponential overshoots
    end
    P_M = P_M ./ sum(P_M);  % Normalize
    
    % Sample the number of firings for each cell (multinomial sampling)
    j2 = mnrnd(NumReacts, P_M);

    % Get state variables modified by this reaction
    cS = find(MRM(:,j1)');
    scS = size(cS);
    SiM2 = SiM;  % Backup current state

    % Apply changes to SiM
    if (scS(2) > 0 && ~any(isnan(j2)))
        SiM(:,cS) = SiM(:,cS) + (MRM(cS,j1) * j2)';

        % Check and correct for any negative values (in case of overreaction)
        while sum(any(SiM(:,cS) < 0)) > 0
            [aN,~] = find(SiM(:,cS) < 0);
            
            % Estimate discrepancy and redistribute reaction events
            diffC = sum(abs(SiM(aN,cS(1))));
            diffCI = floor(diffC);
            diffCd = diffC - diffCI;

            % Attempt compensation between compartments
            if length(cS) > 1
                SiM(aN,cS(2)) = SiM(aN,cS(2)) + SiM(aN,cS(1));
            end
            SiM(aN,cS(1)) = 0;

            % Retry reaction assignment if possible
            if sum(SiM(:,IndS(1,j1))) > 0
                SiM11 = SiM1; SiM11(:,1:NumS) = SiM;
                ReacsS = floor(SiM11(:,IndS(1,j1))) ./ sum(floor(SiM11(:,IndS(1,j1))));
                P_M = 1 - exp(-A1(j1) * dt * ReacsS);
                if sum(P_M) > 0.1
                    P_M = SiM11(:,IndS(1,j1)) ./ sum(SiM11(:,IndS(1,j1)));
                end
                P_M = P_M ./ sum(P_M);
                j3 = mnrnd(diffCI, P_M);
                if ~any(isnan(j3))
                    SiM(:,cS) = SiM(:,cS) + (MRM(cS,j1) * j3)';
                end
            end

            % Add remaining fractional events
            if diffCd > 0 && sum(SiM(:,IndS(1,j1))) > 0
                j3 = mnrnd(1, P_M);
                if ~any(isnan(j3))
                    SiM(:,cS) = SiM(:,cS) + diffCd * (MRM(cS,j1) * j3)';
                end
            end
        end

        % Compute net change from reaction (used for downstream processes)
        j2 = ((SiM(:,cS(1)) - SiM2(:,cS(1))) .* MRM(cS(1),j1))';

        % Update replication cycles if reaction affects replicated genomes
        if j1 ~= Vir_Prod
            SiMC = RepCycles(j2, NumS, maxCycle, A1, j1, IndS, Pos_syn, SiMC, cS, MRM, dt, SiM);
        end

        % Special handling for RNA import into nucleus and addition of new RNPs
        if j1 == RNA_nuc
            if ~any(isnan(j2))
                % Update RNP count globally and per MOI based on imported RNA
                Si(1,nRNP) = Si(1,nRNP) + sum((V0(find(j2),:) > 0) .* j2', 1);
                SiM(:,nRNP) = SiM(:,nRNP) + (V0(find(j2),:) > 0) .* j2';

                IndJ = find(j2); j2 = j2(IndJ);
                for sj = 1:size(j2',1)
                    % Update replication cycle structure for imported RNPs
                    SiMC2 = reshape(SiMC(IndJ(sj),:), NumS, maxCycle)';
                    SiMC2(1,nRNP) = SiMC2(1,nRNP) + (V0(IndJ(sj),:) > 0) .* j2(sj);
                    SiMC(IndJ(sj),:) = reshape(SiMC2', 1, NumS*maxCycle);
                end
            end
        end
    end
end

end
