% This function performs a deterministic update of protein-bound complexes
% following a set of reactions (Rs) involving a specific protein species (Prots_reac).
% It calculates how many complexes dissociate from the protein and updates the
% multivariant state matrix SiM accordingly.

function SiM = ReactProts(Rs,A1,SiM1,SiM,dt,Prots_reac,Params_reac,Si)

    % Compute reaction probabilities based on protein abundance and kinetics
    P_M = 1 - exp(-sum(A1(Rs)) .* dt .* SiM(:,Prots_reac) .* Params_reac ./ sum(SiM(:,Prots_reac) .* Params_reac));

    % If reaction probabilities are too high, use relative abundance instead
    if sum(P_M) > 0.1
        P_M = SiM(:,Prots_reac) .* Params_reac ./ sum(SiM(:,Prots_reac) .* Params_reac);
    end
    P_M = P_M ./ sum(P_M);  % Normalize probabilities

    % Sample how many protein complexes dissociate
    j2 = mnrnd(sum(SiM(:,Prots_reac)) - Si(Prots_reac), P_M);

    % Update SiM if sampling was successful
    if ~any(isnan(j2))
        SiM(:,Prots_reac) = SiM(:,Prots_reac) - j2';
    end

    % Correct potential negative values due to stochastic fluctuations
    count = 0;
    while sum(SiM(:) < 0) > 0
        count = count + 1;
        if count > 100
            break
        end

        % Restore negatives to zero and resample missing amounts
        SiM2 = SiM;
        SiM(SiM < 0) = 0;
        diffP = sum(SiM(:,Prots_reac) - SiM2(:,Prots_reac));
        j2 = mnrnd(diffP, P_M);

        if ~any(isnan(j2))
            SiM(:,Prots_reac) = SiM(:,Prots_reac) - j2';
        end
    end

end
