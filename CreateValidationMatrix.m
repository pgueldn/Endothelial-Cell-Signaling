function Exp_validation=CreateValidationMatrix(speciesNames, numPert)

a=find(strcmp(speciesNames, 'TGFb1'));
if size(a)==1
speciesNames{a}='TGFb';
end

Stress = find(strcmp(speciesNames, 'Stress'));
Wss = find(strcmp(speciesNames, 'Wss'));
Hypoxia = find(strcmp(speciesNames, 'Hypoxia'));
TNFa = find(strcmp(speciesNames, 'TNFa'));
TSP1 = find(strcmp(speciesNames, 'TSP1'));

Proliferation = find(strcmp(speciesNames, 'Proliferation'));
Migration = find(strcmp(speciesNames, 'Migration'));
Apoptosis = find(strcmp(speciesNames, 'Apoptosis'));
EndMT = find(strcmp(speciesNames, 'EndMT'));
MCP1 = find(strcmp(speciesNames, 'MCP1'));
HIF1a = find(strcmp(speciesNames, 'HIF1a'));
NFkB = find(strcmp(speciesNames, 'NFkB'));
VEGF = find(strcmp(speciesNames, 'VEGF'));
VEGFR = find(strcmp(speciesNames, 'VEGFR'));
IL6 = find(strcmp(speciesNames, 'IL6'));
ppET1 = find(strcmp(speciesNames, 'ppET1'));
FAK = find(strcmp(speciesNames, 'FAK'));
pSmad23 = find(strcmp(speciesNames, 'pSmad23'));
pSmad15 = find(strcmp(speciesNames, 'pSmad15'));
p38 = find(strcmp(speciesNames, 'p38'));
ERK = find(strcmp(speciesNames, 'ERK'));
JNK = find(strcmp(speciesNames, 'JNK'));
mTOR = find(strcmp(speciesNames, 'mTOR'));
PKC = find(strcmp(speciesNames, 'PKC'));
DAG = find(strcmp(speciesNames, 'DAG'));
RhoA = find(strcmp(speciesNames, 'RhoA'));
ROCK = find(strcmp(speciesNames, 'ROCK'));
PI3K = find(strcmp(speciesNames, 'PI3K'));
Akt = find(strcmp(speciesNames, 'Akt'));
eNOS = find(strcmp(speciesNames, 'eNOS'));
ET1 = find(strcmp(speciesNames, 'ET1'));
TGFb = find(strcmp(speciesNames, 'TGFb'));
MLCK = find(strcmp(speciesNames, 'MLCK'));
p70S6K = find(strcmp(speciesNames, 'p70S6K'));
Actin = find(strcmp(speciesNames, 'Actin'));
TSP1 = find(strcmp(speciesNames, 'TSP1'));

Exp_validation=0.5*ones(numPert,length(speciesNames)); %a non-zero value so 'no change' can be recorded as 0.

Exp_validation(Stress,[Proliferation, Migration, EndMT, MCP1, NFkB, VEGFR, IL6, TGFb, ppET1, FAK, p38, ERK, JNK, RhoA, PI3K, Akt, eNOS, MLCK, p70S6K, Actin])=1; %Stress increases
Exp_validation(Stress,[])=0; %Stress doesn't change
Exp_validation(Stress,[Apoptosis])=-1; %Stress decreases
Exp_validation(Stress,[])=-0.5; %mixed

Exp_validation(Wss,[Migration, pSmad15, VEGF, VEGFR, FAK, MLCK, ROCK, ERK, JNK, PKC, p38, PI3K, Akt, eNOS, RhoA, NFkB])=1;
Exp_validation(Wss,[])=0;
Exp_validation(Wss,[Proliferation, EndMT, Apoptosis, ppET1, ET1, pSmad23])=-1;
Exp_validation(Wss,[MCP1, TGFb])=-0.5; %mixed

Exp_validation(Hypoxia,[Proliferation, Migration, Apoptosis, EndMT, MCP1, VEGF, IL6, ppET1, ET1, NFkB, TSP1, TGFb, pSmad23,...
    mTOR, Akt, p70S6K, ROCK, p38, JNK, ERK])=1;
Exp_validation(Hypoxia,[])=0;
Exp_validation(Hypoxia,[])=-1;
Exp_validation(Hypoxia,[eNOS])=-0.5; %mixed 

Exp_validation(TNFa-1,[Apoptosis, EndMT, NFkB, IL6, p38, JNK, ERK, ppET1, ET1, RhoA, ROCK, MLCK, MCP1, PI3K, Akt])=1;
Exp_validation(TNFa-1,[])=0; 
Exp_validation(TNFa-1,[Proliferation])=-1;
Exp_validation(TNFa-1,[eNOS, VEGFR])=-0.5; %mixed 

end
