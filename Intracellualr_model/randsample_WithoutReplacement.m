function C = randsample_WithoutReplacement(V,k,P)
    
    P(P < 0) = 0;
    P = P./sum(P);

    if length(V) == k

        C = V;

    else

        C = zeros(1,k);
        
        for I = 1:k
            
            p = find(mnrnd(1,P));
            C(I) = V(p);
            P(p) = []; P = P./sum(P);
            V(p) = [];

        end

    end
    
    C = sort(C);
end
