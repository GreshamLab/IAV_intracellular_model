% This function performs a deterministic update of protein-bound complexes
% across multiple genotypes without explicit reaction kinetics.
% It redistributes the difference between the current and target total counts
% of a given protein species (Prots_reac) across genotypes according to their
% relative abundance, and updates the multivariant state matrix SiM.

function SiM = ReactProtsGenDet(SiM, Prots_reac, Params_reac, Si)

    % Compute the relative contribution of each genotype based on protein abundance and kinetic parameters
    P_M = SiM(:,Prots_reac) .* Params_reac ./ sum(SiM(:,Prots_reac) .* Params_reac);

    % Normalize probabilities across genotypes
    P_M = P_M ./ sum(P_M, 1);

    % Deterministically subtract the excess molecules from each genotype
    % according to their normalized contributions
    SiM(:,Prots_reac) = SiM(:,Prots_reac) - ...
        (sum(SiM(:,Prots_reac)) - Si(Prots_reac)) .* P_M;

    % If the target count is zero, clear all entries for this protein
    if Si(Prots_reac) == 0
        SiM(:,Prots_reac) = 0;
    end

end
