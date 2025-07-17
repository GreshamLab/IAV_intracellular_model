% This function performs a deterministic update of protein-bound complexes
% following a set of reactions (Rs) involving a specific protein species (Prots_reac).
% It calculates the fraction of complexes that dissociate during a time step (dt)
% and updates the multivariant state matrix SiM accordingly, without stochastic sampling.

function SiM = ReactProtsDet(Rs,A1,SiM,dt,Prots_reac,Params_reac,Si)

    % Compute reaction probabilities based on reaction rates and protein abundances
    P_M = 1 - exp(-sum(A1(Rs)) .* dt .* SiM(:,Prots_reac) .* Params_reac ./ ...
                  sum(SiM(:,Prots_reac) .* Params_reac));
    
    % Total probability of dissociation
    P_T = sum(P_M, 1);
    
    % If total probability exceeds 1 (numerical instability), use relative contributions instead
    if P_T > 1
        P_M = SiM(:,Prots_reac) .* Params_reac ./ sum(SiM(:,Prots_reac) .* Params_reac);
    end
    
    % Normalize probabilities across all genotypes
    P_M = P_M ./ sum(P_M, 1);
    
    % Deterministically redistribute the reduction in protein complexes
    % according to the normalized probabilities
    SiM(:,Prots_reac) = SiM(:,Prots_reac) - ...
        (sum(SiM(:,Prots_reac)) - Si(Prots_reac)) .* P_M;

    % If no target complexes remain, zero out this protein's entries
    if Si(Prots_reac) == 0
        SiM(:,Prots_reac) = 0;
    end

end
