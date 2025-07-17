function[Si,SiM,SiMC] = detEulerMOI(A1,MRM,dt,SiM1,SiM,IndS,NumS,maxCycle,Pos_syn,SiMC,MOI,Si,RNA_nuc,V0,nRNP,Rdet)
    %%
    dM1 = (A1(Rdet).*MRM(:,Rdet))';

    P_M = 1 - exp(-A1.*dt.*SiM1(:,IndS(1,:))./sum(SiM1(:,IndS(1,:))));
    P_T = sum(P_M,1);
    P_M(:,P_T > 1) = SiM1(:,IndS(1,P_T > 1))./sum(SiM1(:,IndS(1,P_T > 1)));
    P_M = P_M./sum(P_M,1); P_M(isnan(P_M)) = 0;
    dSiM = P_M(:,Rdet)*dM1.*dt;
    SiM = SiM + dSiM;

    if ismember(RNA_nuc,Rdet)
        SiM(:,nRNP) = SiM(:,nRNP) + (V0>0).*P_M(:,nRNP).*A1(RNA_nuc).*dt;
        Si(1,nRNP) = Si(1,nRNP) + sum((V0>0).*P_M(:,nRNP).*A1(RNA_nuc).*dt,1);
    end
    
    [Mneg,Sneg] = find(SiM<0);

    for sN = 1:length(Sneg)

        dSiM(:,Sneg(sN)) = -SiM(Mneg(sN),Sneg(sN));
        SiM(Mneg(sN),Sneg(sN)) = 0;
        SiM(:,Sneg(sN)) = SiM(:,Sneg(sN)) - dSiM(:,Sneg(sN)).*SiM(:,Sneg(sN))./sum(SiM(:,Sneg(sN)),1);

    end
%%
    for m = 1:MOI
        
        SiMC = detEulerCycles(A1,MRM,P_M,NumS,maxCycle,IndS,Pos_syn,SiMC,dt,m,Rdet,SiM);
        
        if ismember(RNA_nuc,Rdet)
            SiM(m,nRNP) = SiM(m,nRNP) + (V0(m,:)>0).*P_M(m,nRNP).*A1(RNA_nuc).*dt;
            SiMC2 = reshape(SiMC(m,:),NumS,maxCycle)';
            SiMC2(1,nRNP) = SiMC2(1,nRNP) + (V0(m,:)>0).*P_M(m,nRNP).*A1(RNA_nuc).*dt;
            SiMC(m,:) = reshape(SiMC2(:,1:NumS)',1,NumS*maxCycle);
        end

    end
  %%  
end
