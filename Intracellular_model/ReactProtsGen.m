% This function performs a generalized update of protein-bound complexes.
% It distributes a specified number of protein complexes (Prots_reac)
% according to a normalized weight vector based on binding parameters.
% This is typically used for generic regulatory protein interactions.

function SiM = ReactProtsGen(SiM1, SiM, Prots_reac, Params_reac, Si)

    % Compute distribution probabilities based on protein abundance and interaction strength
    P_M = SiM(:,Prots_reac) .* Params_reac ./ sum(SiM(:,Prots_reac) .* Params_reac);
    P_M = P_M ./ sum(P_M);  % Normalize to form a probability distribution

    % Sample how many complexes are released from proteins
    j2 = mnrnd(sum(SiM(:,Prots_reac)) - Si(Prots_reac), P_M);

    % Apply the sampled changes to the system if no NaNs were produced
    if ~any(isnan(j2))
        SiM(:,Prots_reac) = SiM(:,Prots_reac) - j2';
    end

    % Ensure no negative values remain after update (numerical robustness)
    count = 0;
    while sum(SiM(:) < 0) > 0
        count = count + 1;
        if count > 100
            break
        end

        % Reset negative values and re-distribute the discrepancy
        SiM2 = SiM;
        SiM(SiM < 0) = 0;

        diffP = sum(SiM(:,Prots_reac) - SiM2(:,Prots_reac));
        j2 = mnrnd(diffP, P_M);

        if ~any(isnan(j2))
            SiM(:,Prots_reac) = SiM(:,Prots_reac) - j2';
        end
    end

end
