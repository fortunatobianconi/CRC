% ODEs of the p38 MAPK signaling pathway taken from:
% Peng, Huiming, et al. 
% "Characterization of p38 MAPK isoforms for drug resistance study using systems biology approach." 
% Bioinformatics 30.13 (2014): 1899-1907.

function dy=MM_model(t,y,parameter)

k=parameter;

dy = zeros(40,1);

%%----------------%%
%    Parameters    %
%%----------------%%

%---------Module 1-------------------
k_GFR=k(1);
k_pGFR=k(2);
k_Shc_pGFR=k(3);
k_pShc=k(4);
%---------Module 2-------------------
k_RasGDP_pShc=k(5);
k_RasGTP_pERK12=k(6);
%---------Module 3-------------------
k_IRS1_pGFR=k(7);
k_pIRS1_pERK12=k(8);
k_pIRS1_pAKT=k(9);
k_PI3K_pIRS1=k(10);
k_pPI3K=k(11);
%---------Module 4-------------------
k_MKK47_RasGTP=k(12);
k_MKK47_pPI3K=k(13);
k_pMKK47=k(14);
k_p38_pMKK47=k(15);
k_pp38=k(16);
%---------Module 5-------------------
k_PDK1_pPI3K=k(17); 
k_pPDK1=k(18); 
%---------Module 6-------------------
k_AKT_pPDK1=k(19); 
k_AKT_pPI3K=k(20);
k_AKT_RasGTP=k(21);
k_pAKT=k(22);
%---------Module 7------------------
k_mTOR_pAKT=k(23); 
k_pmTOR=k(24); 
%---------Module 8-------------------
k_Raf1_RasGTP=k(25);
k_Raf1_pPI3K=k(26); 
k_pRaf1_pAKT=k(27);
%---------Module 9-------------------
k_MEK12_pRaf1=k(28); 
k_pMEK12=k(29);
%---------Module 10-------------------
k_ERK12_pMEK12=k(30); 
k_ERK12_pp38=k(31); 
k_pERK12=k(32);
%---------Module 11------------------
k_P70S6K_pERK12=k(33); 
k_P70S6K_pmTOR=k(34);
k_pP70S6K=k(35);
%---------Module 12------------------
k_JNK_pMKK47=k(36);
k_pJNK=k(37);
k_cJUN_pJNK=k(38);
k_pcJUN=k(39);
%---------Module 13------------------
k_BCLXL_pp38=k(40);
k_BCLXL_pJNK=k(41);
%---------Module 14------------------
k_BAX_pJNK=k(42);
k_BAX_pp38=k(43);
k_BAX=k(44);
%---------Module 15------------------
k_IKK_pMKK47=k(45);
k_IKK_pAKT=k(46);
k_pIKK=k(47);
k_NFkB_pIKK=k(48);
k_pNFkB=k(49);
%---------Module 16------------------
k_PARP_pcJUN=k(50);
k_PARP_BAX=k(51);
k_cPARP_BCLXL=k(52);
k_cPARP_pNFkB=k(53);


%%------------------------------%%
%   Concentration of Proteins    %
%%------------------------------%%

%---------Module 1------------------- 
GFR=y(1);
pGFR=y(2);
Shc=y(3);
pShc=y(4);
%---------Module 2-------------------
RasGDP=y(5);
RasGTP=y(6);
%---------Module 3-------------------
IRS1=y(7);
pIRS1=y(8);
PI3K=y(9);
pPI3K=y(10);
%---------Module 4-------------------
MKK47=y(11);
pMKK47=y(12);
p38=y(13);
pp38=y(14);
%---------Module 5-------------------
PDK1=y(15);
pPDK1=y(16);
%---------Module 6-------------------
AKT=y(17);
pAKT=y(18);
%---------Module 7-------------------
mTOR=y(19);
pmTOR=y(20);
%---------Module 8-------------------
Raf1=y(21);
pRaf1=y(22);
%---------Module 9-------------------
MEK12=y(23);
pMEK12=y(24);
%---------Module 10-------------------
ERK12=y(25);
pERK12=y(26);
%---------Module 11-------------------
P70S6K=y(27);
pP70S6K=y(28);
%---------Module 12-------------------
JNK=y(29);
pJNK=y(30);
cJUN=y(31);
pcJUN=y(32);
%---------Module 13------------------
BCLXL=y(33);
%---------Module 14------------------
BAX=y(34);
%---------Module 15-------------------
IKK=y(35);
pIKK=y(36);
NFkB=y(37);
pNFkB=y(38);
%---------Module 16-------------------
PARP=y(39);
cPARP=y(40);


%%-----------------------------------%%
%              Equations              %
%%-----------------------------------%%

%---------Module 1-------------------
d_GFR=-k_GFR*GFR+k_pGFR*pGFR;
d_pGFR=k_GFR*GFR-k_pGFR*pGFR;
d_Shc=-k_Shc_pGFR*Shc*pGFR+k_pShc*pShc;
d_pShc=k_Shc_pGFR*Shc*pGFR-k_pShc*pShc;
%---------Module 2-------------------
d_RasGDP=-k_RasGDP_pShc*RasGDP*pShc+k_RasGTP_pERK12*RasGTP*pERK12;
d_RasGTP= k_RasGDP_pShc*RasGDP*pShc-k_RasGTP_pERK12*RasGTP*pERK12;
%---------Module 3-------------------
d_IRS1=-k_IRS1_pGFR*IRS1*pGFR+k_pIRS1_pERK12*pIRS1*pERK12+k_pIRS1_pAKT*pIRS1*pAKT;
d_pIRS1=k_IRS1_pGFR*IRS1*pGFR-k_pIRS1_pERK12*pIRS1*pERK12-k_pIRS1_pAKT*pIRS1*pAKT;
d_PI3K=-k_PI3K_pIRS1*PI3K*pIRS1+k_pPI3K*pPI3K;
d_pPI3K=k_PI3K_pIRS1*PI3K*pIRS1-k_pPI3K*pPI3K;
%---------Module 4-------------------
d_MKK47=-k_MKK47_RasGTP*MKK47*RasGTP-k_MKK47_pPI3K*MKK47*pPI3K+k_pMKK47*pMKK47;
d_pMKK47=k_MKK47_RasGTP*MKK47*RasGTP+k_MKK47_pPI3K*MKK47*pPI3K-k_pMKK47*pMKK47;
d_p38=-k_p38_pMKK47*p38*pMKK47+k_pp38*pp38;
d_pp38=k_p38_pMKK47*p38*pMKK47-k_pp38*pp38;
%---------Module 5-------------------
d_PDK1=-k_PDK1_pPI3K*PDK1*pPI3K+k_pPDK1*pPDK1;
d_pPDK1=k_PDK1_pPI3K*PDK1*pPI3K-k_pPDK1*pPDK1;
%---------Module 6-------------------
d_AKT=-k_AKT_pPDK1*AKT*pPDK1-k_AKT_pPI3K*AKT*pPI3K-k_AKT_RasGTP*AKT*RasGTP+k_pAKT*pAKT;
d_pAKT=k_AKT_pPDK1*AKT*pPDK1+k_AKT_pPI3K*AKT*pPI3K+k_AKT_RasGTP*AKT*RasGTP-k_pAKT*pAKT;
%---------Module 7-------------------
d_mTOR=-k_mTOR_pAKT*mTOR*pAKT+k_pmTOR*pmTOR;
d_pmTOR=k_mTOR_pAKT*mTOR*pAKT-k_pmTOR*pmTOR;
%---------Module 8-------------------
d_Raf1=-k_Raf1_RasGTP*Raf1*RasGTP-k_Raf1_pPI3K*Raf1*pPI3K+k_pRaf1_pAKT*pRaf1*pAKT;
d_pRaf1=k_Raf1_RasGTP*Raf1*RasGTP+k_Raf1_pPI3K*Raf1*pPI3K-k_pRaf1_pAKT*pRaf1*pAKT;
%---------Module 9-------------------
d_MEK12=-k_MEK12_pRaf1*MEK12*pRaf1+k_pMEK12*pMEK12;
d_pMEK12=k_MEK12_pRaf1*MEK12*pRaf1-k_pMEK12*pMEK12;
%---------Module 10-------------------
d_ERK12=-k_ERK12_pMEK12*ERK12*pMEK12-k_ERK12_pp38*ERK12*pp38+k_pERK12*pERK12;
d_pERK12=k_ERK12_pMEK12*ERK12*pMEK12+k_ERK12_pp38*ERK12*pp38-k_pERK12*pERK12;
%---------Module 11-------------------
d_P70S6K=-k_P70S6K_pERK12*P70S6K*pERK12-k_P70S6K_pmTOR*P70S6K*pmTOR+k_pP70S6K*pP70S6K;
d_pP70S6K=k_P70S6K_pERK12*P70S6K*pERK12+k_P70S6K_pmTOR*P70S6K*pmTOR-k_pP70S6K*pP70S6K;
%---------Module 12-------------------
d_JNK=-k_JNK_pMKK47*JNK*pMKK47+k_pJNK*pJNK;
d_pJNK=k_JNK_pMKK47*JNK*pMKK47-k_pJNK*pJNK;
d_cJUN=-k_cJUN_pJNK*cJUN*pJNK+k_pcJUN*pcJUN;
d_pcJUN=k_cJUN_pJNK*cJUN*pJNK-k_pcJUN*pcJUN;
%---------Module 13------------------
d_BCLXL=k_BCLXL_pp38*pp38-k_BCLXL_pJNK*BCLXL*pJNK;
%---------Module 14------------------
d_BAX=k_BAX_pJNK*pJNK+k_BAX_pp38*pJNK*pp38-k_BAX*BAX;
%---------Module 15-------------------
d_IKK=-k_IKK_pMKK47*IKK*pMKK47-k_IKK_pAKT*IKK*pAKT+k_pIKK*pIKK;
d_pIKK=k_IKK_pMKK47*IKK*pMKK47+k_IKK_pAKT*IKK*pAKT-k_pIKK*pIKK;
d_NFkB=-k_NFkB_pIKK*NFkB*pIKK+k_pNFkB*pNFkB;
d_pNFkB=k_NFkB_pIKK*NFkB*pIKK-k_pNFkB*pNFkB;
%---------Module 16-------------------
d_PARP=-k_PARP_pcJUN*PARP*pcJUN-k_PARP_BAX*PARP*BAX+k_cPARP_BCLXL*cPARP*BCLXL+k_cPARP_pNFkB*cPARP*pNFkB;
d_cPARP=k_PARP_pcJUN*PARP*pcJUN+k_PARP_BAX*PARP*BAX-k_cPARP_BCLXL*cPARP*BCLXL-k_cPARP_pNFkB*cPARP*pNFkB;



%%-----------------------------------%%
%           Differential Terms        %
%%-----------------------------------%%

%---------Module 1-------------------
dy(1)=d_GFR;
dy(2)=d_pGFR;
dy(3)=d_Shc;
dy(4)=d_pShc;
%---------Module 2-------------------
dy(5)=d_RasGDP;
dy(6)=d_RasGTP;
%---------Module 3-------------------
dy(7)=d_IRS1;
dy(8)=d_pIRS1;
dy(9)=d_PI3K;
dy(10)=d_pPI3K;
%---------Module 4-------------------
dy(11)=d_MKK47;
dy(12)=d_pMKK47;
dy(13)=d_p38;
dy(14)=d_pp38;
%---------Module 5-------------------
dy(15)=d_PDK1;
dy(16)=d_pPDK1;
%---------Module 6-------------------
dy(17)=d_AKT;
dy(18)=d_pAKT;
%---------Module 7-------------------
dy(19)=d_mTOR;
dy(20)=d_pmTOR;
%---------Module 8-------------------
dy(21)=d_Raf1;
dy(22)=d_pRaf1;
%---------Module 9-------------------
dy(23)=d_MEK12;
dy(24)=d_pMEK12;
%---------Module 10-------------------
dy(25)=d_ERK12;
dy(26)=d_pERK12;
%---------Module 11-------------------
dy(27)=d_P70S6K;
dy(28)=d_pP70S6K;
%---------Module 12-------------------
dy(29)=d_JNK;
dy(30)=d_pJNK;
dy(31)=d_cJUN;
dy(32)=d_pcJUN;
%---------Module 13------------------
dy(33)=d_BCLXL;
%---------Module 14------------------
dy(34)=d_BAX;
%---------Module 15-------------------
dy(35)=d_IKK;
dy(36)=d_pIKK;
dy(37)=d_NFkB;
dy(38)=d_pNFkB;
%---------Module 16-------------------
dy(39)=d_PARP;
dy(40)=d_cPARP;