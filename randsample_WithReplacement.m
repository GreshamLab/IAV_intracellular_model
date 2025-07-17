% This function performs random sampling without replacement from a weighted distribution.
% It draws `k` unique elements from the vector `V`, using the corresponding probabilities in `P`.
% The sampled elements are returned in sorted order.

function C = randsample_WithReplacement(V, k, P)

    % Sanitize input probabilities: remove negatives and normalize
    P(P < 0) = 0;
    P = P ./ sum(P);

    % If all elements must be selected, return the full set
    if length(V) == k
        C = V;

    else
        % Initialize output vector
        C = zeros(1, k);

        % Iteratively sample elements without replacement
        for I = 1:k
            % Draw one index based on multinomial sampling
            p = find(mnrnd(1, P));

            % Store sampled element
            C(I) = V(p);

            % Remove the selected element and re-normalize probabilities
            P(p) = [];
            P = P ./ sum(P);
            V(p) = [];
        end
    end

    % Return sampled elements in sorted order
    C = sort(C);
end
