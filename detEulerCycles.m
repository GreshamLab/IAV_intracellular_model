function SiMC = detEulerCycles(A1,MRM,P_M,NumS,maxCycle,IndS,Pos_syn,SiMC,dt,m,Rdet,SiM)
%%
dM1 = (A1(Rdet).*MRM(:,Rdet))'.*P_M(m,Rdet)';

SiMC2 = ones(maxCycle,NumS+1);
SiMC2(:,1:NumS) = reshape(SiMC(m,:),NumS,maxCycle)';
P_C = 1 - exp(-A1.*dt.*SiMC2(:,IndS(1,:))./sum(SiMC2(:,IndS(1,:))));
P_T = sum(P_C,1);
P_C(:,P_T > 1) = SiMC2(:,IndS(1,P_T > 1))./sum(SiMC2(:,IndS(1,P_T > 1))); 
P_C = P_C./sum(P_C,1); P_C(isnan(P_C)) = 0;

P_C(2:end,Pos_syn) = P_C(1:end-1,Pos_syn); P_C(1,Pos_syn) = 0;
dSiMC = P_C(:,Rdet)*dM1.*dt;
SiMC2(:,1:NumS) = SiMC2(:,1:NumS) + dSiMC;

[Mneg,Sneg] = find(SiMC2<0);

for sN = 1:length(Sneg)
    negs = find(SiMC2(:,Sneg(sN)) < 0);
    SiMC2(negs,Sneg(sN)) = 0;
    realSM = SiM(m,Sneg(sN));
    if sum(SiMC2(:,Sneg(sN))) > 0
        SiMC2(:,Sneg(sN)) = realSM.*(SiMC2(:,Sneg(sN))./sum(SiMC2(:,Sneg(sN))));
    end
    
end

% if sum(abs(SiM(m,1:35) - sum(SiMC2(:,1:35),1)) > 1) > 0
%     SiM(m,1:35) - sum(SiMC2(:,1:35),1)
%     m
%     brea; end
% 
% if sum(SiMC2(:)<0) > 0; brea; end

SiMC(m,:) = reshape(SiMC2(:,1:NumS)',1,NumS*maxCycle);
end