function [V0,Infecting_virions_locusC,Infecting_virions_mutsC,Infecting_virions_numMC] = Inf_virions(NumG,Length_muts,Cell,Initial_virions_segments,Seg_coor,Infecting_virions_locus,Infecting_virions_muts,Infecting_virions_numM,a,Initial_virions_cycle,IAV_segments)

Infecting_virions_numMC = zeros(1,NumG);
Infecting_virions_locusC = zeros(1,Length_muts);
Infecting_virions_mutsC = zeros(1,Length_muts);

virions = sum(Initial_virions_segments(:,:,Cell),2); virions = length(find(virions>0));
if virions(1) > 0
    virion1 = randi(virions(1),1,1);
else
    virion1 = 1;
end

V0 = Initial_virions_segments(virion1,:,Cell);

% Infecting_virions_numMC(c,:) = Infecting_virions_numM(V0(c,:)==(1:maxMOI),:,Cells);
c1 = 0;
fS = find(V0);
%%
for j = 1:length(fS)

    I = fS(j);
    Infecting_virions_numMC(1,I) = Infecting_virions_numM(V0(1,I),I,Cell);

    if Infecting_virions_numMC(1,I) == 0; continue; end
    ic = c1 + 1;
    D = zeros(Length_muts,2);
    D(:,1) = (Infecting_virions_locus(V0(1,I),:,Cell) - Seg_coor(1,I)+1) > 0;
    D(:,2) = (Seg_coor(2,I)+1 - Infecting_virions_locus(V0(1,I),:,Cell)) > 0;
    
    M1 = find(prod(D,2));
    c1 = c1 + length(M1);
    Infecting_virions_locusC(1,ic:c1) = Infecting_virions_locus(V0(1,I),M1,Cell);
    Infecting_virions_mutsC(1,ic:c1) = Infecting_virions_muts(V0(1,I),M1,Cell);

end
%%
Muts = poissrnd(a.*IAV_segments.*Initial_virions_cycle(virion1,:,Cell));
fMuts = find(Muts); numM = length(fMuts);

for s = 1:numM
    
    seg = fMuts(s);
    locus = randi([Seg_coor(1,seg),Seg_coor(2,seg)],1,Muts(seg));
    mut = randi([2,4],1,Muts(seg));
    
    iniC = nansum(Infecting_virions_numMC(1,:)) + 1;
    finC = nansum(Infecting_virions_numMC(1,:)) + Muts(seg);
    Infecting_virions_locusC(1,iniC:finC) = locus;
    Infecting_virions_mutsC(1,iniC:finC) = mut;
    Infecting_virions_numMC(1,seg) = Infecting_virions_numMC(1,seg) + Muts(seg);
    
end


end