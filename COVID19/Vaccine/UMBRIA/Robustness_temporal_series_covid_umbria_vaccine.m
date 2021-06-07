tic;

init_umbria_vaccine

% Load parameters (set to 1 or mode previous step)
%load('mode_parameters_2step_LINLOG_Nr10_Ns100000_UMBRIA_I2I3DE_multivariate_1SET_2may.mat');
pN=ones(1,31);


% Generate perturbation matrix usign latin hypercube sampling
Nsample=100000; % number of values for each parameters
Nr=10; % number of samples

Results=cell(Nr,4);

poolobj=parpool(4);
for k=1:Nr
    k
    
    %uncomment from second step onwards
    %pN=mode_parameters_2step_LINLOG(k,:);
 
    %bounds of transmission rate parameters 
    LBpi_be1=0.01/pN(1);   
    UBpi_be1=1/pN(1);     
 
    LBpi_b01=0.01/pN(2);  
    UBpi_b01=1/pN(2);    
    
    LBpi_b1=0.001/pN(3);    
    UBpi_b1=1/pN(3);    
 
    LBpi_b2=0.001/pN(4);    
    UBpi_b2=1/pN(4);  
 
    LBpi_b3=0.001/pN(5);   
    UBpi_b3=1/pN(5);     
    
    %bounds for fractions FracSev
    LBpi_FracSev=0.01/pN(6);  
    UBpi_FracSev=0.08/pN(6);   
 
    %bounds for fractions FracCrit
    LBpi_FracCrit=0.001/pN(7);   
    UBpi_FracCrit=0.02/pN(7);    
 
    %bounds for fractions FracAsym
    LBpi_FracAsym = 0.2/pN(8);   
    UBpi_FracAsym = 0.7/pN(8);                  
    
    %bounds for incubation period, mild infection period and asymptomatic
    %infection period (IncubPeriod, DurMildInf,DurAsym)
    LBpi_Incub=4/pN(9);    
    UBpi_Incub=6/pN(9);   
 
    LBpi_DurMild=5/pN(10);  
    UBpi_DurMild=30/pN(10);  
 
    LBpi_DurAsym=5/pN(11);   
    UBpi_DurAsym=20/pN(11);             
 
    %bounds for duration severe infection (DurHosp)
    LBpi_DurHosp=4/pN(12);   
    UBpi_DurHosp=30/pN(12);          
 
    %bounds for ICU stay (TimeICUDeath)
    LBpi_ICU=4/pN(13);   
    UBpi_ICU=30/pN(13);       
 
    %bounds for ProbDeath
    LBpi_ProbDeath=10/pN(14);   
    UBpi_ProbDeath=90/pN(14);  
 
    %bounds for PresymPeriod
    LBpi_Presym=0.5/pN(15);  
    UBpi_Presym=0.9/pN(15);            
    
    %bounds for s0_0
    LBpi_s0_0=0.5/pN(16);  
    UBpi_s0_0=1.5/pN(16);  
    
    %bounds for s0_1
    LBpi_s0_1=0.1/pN(17);  
    UBpi_s0_1=0.9/pN(17); 
 
 
    %bounds for s0_2
    LBpi_s0_2=0.1/pN(18);  
    UBpi_s0_2=0.9/pN(18); 
    
     
    %bounds for s0_3
    LBpi_s0_3=0.5/pN(19);  
    UBpi_s0_3=1.5/pN(19);  
     
    
    %bounds for s0_4
    LBpi_s0_4=0.5/pN(20);  
    UBpi_s0_4=1.5/pN(20);
 
    %bounds for s0_5
    LBpi_s0_5=0.1/pN(21);  
    UBpi_s0_5=0.9/pN(21);    
    
    %bounds for s0_6
    LBpi_s0_6=0.1/pN(22);  
    UBpi_s0_6=0.9/pN(22);     
    
    %bounds for fraction FracCrit2
    LBpi_FracCrit2=0.002/pN(23);   
    UBpi_FracCrit2=0.02/pN(23);
    
    %bounds for fraction FracCrit3
    LBpi_FracCrit3=0.002/pN(24);   
    UBpi_FracCrit3=0.02/pN(24);
    
    %bounds for fraction FracCrit4
    LBpi_FracCrit4=0.002/pN(25);   
    UBpi_FracCrit4=0.02/pN(25);
    
    %bounds for fraction FracCrit5
    LBpi_FracCrit5=0.002/pN(26);   
    UBpi_FracCrit5=0.02/pN(26);
    
    %bounds for be2
    LBpi_be2=0.01/pN(27);   
    UBpi_be2=1/pN(27);
    
    %bounds for b02
    LBpi_b02=0.01/pN(28);   
    UBpi_b02=1/pN(28);    
  
    %bounds for probability of death of hospitalized patients
    LBpi_ProbDeathH=10/pN(29);   
    UBpi_ProbDeathH=90/pN(29);
    
    %bounds for n
    LBpi_nHILL=1/pN(30);   
    UBpi_nHILL=100/pN(30);
    
    %bounds for K
    LBpi_KHILL=1/pN(31);   
    UBpi_KHILL=10^5/pN(31);

    %perturbation array for transmission rate parameters (be,b0,b1,b2,b3)
    lb_be1=log(LBpi_be1*ones(1,1));
    ub_be1=log(UBpi_be1*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_be1-lb_be1),lb_be1);
    Perturbation_be1(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_be=LBpi_be+(UBpi_be-LBpi_be).*lhs;
    
    
    lb_b01=log(LBpi_b01*ones(1,1));
    ub_b01=log(UBpi_b01*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_b01-lb_b01),lb_b01);
    Perturbation_b01(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_b0=LBpi_b0+(UBpi_b0-LBpi_b0).*lhs;
    
   
    lb_b1=log(LBpi_b1*ones(1,1));
    ub_b1=log(UBpi_b1*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_b1-lb_b1),lb_b1);
    Perturbation_b1(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_b1=LBpi_b1+(UBpi_b1-LBpi_b1).*lhs;
    
    lb_b2=log(LBpi_b2*ones(1,1));
    ub_b2=log(UBpi_b2*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_b2-lb_b2),lb_b2);
    Perturbation_b2(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_b2=LBpi_b2+(UBpi_b2-LBpi_b2).*lhs;
    
    
    lb_b3=log(LBpi_b3*ones(1,1));
    ub_b3=log(UBpi_b3*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_b3-lb_b3),lb_b3);
    Perturbation_b3(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_b3=LBpi_b3+(UBpi_b3-LBpi_b3).*lhs;
    
    %perturbation array for fraction parameters (FracSevere, FracCritical, FracAsym)
    lb_FracSev=log(LBpi_FracSev*ones(1,1));
    ub_FracSev=log(UBpi_FracSev*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_FracSev-lb_FracSev),lb_FracSev);
    Perturbation_FracSev(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_FracSev=LBpi_FracSev+(UBpi_FracSev-LBpi_FracSev).*lhs;
    
    
    lb_FracCrit=log(LBpi_FracCrit*ones(1,1));
    ub_FracCrit=log(UBpi_FracCrit*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_FracCrit-lb_FracCrit),lb_FracCrit);
    Perturbation_FracCrit(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_FracCrit=LBpi_FracCrit+(UBpi_FracCrit-LBpi_FracCrit).*lhs;
    

    %lb_FracAsym=log(LBpi_FracAsym*ones(1,1));
    %ub_FracAsym=log(UBpi_FracAsym*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    %PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_FracAsym-lb_FracAsym),lb_FracAsym);
    %Perturbation_FracAsym(:,:)=exp(PerturbationLog(:,:));
    Perturbation_FracAsym=LBpi_FracAsym+(UBpi_FracAsym-LBpi_FracAsym).*lhs;

        
    %perturbation array for incubation, mild infection and asymptomatic
    %infection period (IncubPeriod, DurMildInf,DurAsym)
    %lb_Incub=log(LBpi_Incub*ones(1,1));
    %ub_Incub=log(UBpi_Incub*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    %PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_Incub-lb_Incub),lb_Incub);
    %Perturbation_Incub(:,:)=exp(PerturbationLog(:,:));
    Perturbation_Incub=LBpi_Incub+(UBpi_Incub-LBpi_Incub).*lhs;

    
    %lb_DurMild=log(LBpi_DurMild*ones(1,1));
    %ub_DurMild=log(UBpi_DurMild*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    %PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_DurMild-lb_DurMild),lb_DurMild);
    %Perturbation_DurMild(:,:)=exp(PerturbationLog(:,:));
    Perturbation_DurMild=LBpi_DurMild+(UBpi_DurMild-LBpi_DurMild).*lhs;
    
    %lb_DurAsym=log(LBpi_DurAsym*ones(1,1));
    %ub_DurAsym=log(UBpi_DurAsym*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    %PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_DurAsym-lb_DurAsym),lb_DurAsym);
    %Perturbation_DurAsym(:,:)=exp(PerturbationLog(:,:));
    Perturbation_DurAsym=LBpi_DurAsym+(UBpi_DurAsym-LBpi_DurAsym).*lhs;
    

    %perturbation array duration of severe infection (DurHosp)
    %lb_DurHosp=log(LBpi_DurHosp*ones(1,1));
    %ub_DurHosp=log(UBpi_DurHosp*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    %PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_DurHosp-lb_DurHosp),lb_DurHosp);
    %Perturbation_DurHosp(:,:)=exp(PerturbationLog(:,:));
    Perturbation_DurHosp=LBpi_DurHosp+(UBpi_DurHosp-LBpi_DurHosp).*lhs;

    
    %perturbation array ICU stay (TimeICUDeath)
    %lb_ICU=log(LBpi_ICU*ones(1,1));
    %ub_ICU=log(UBpi_ICU*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    %PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_ICU-lb_ICU),lb_ICU);
    %Perturbation_ICU(:,:)=exp(PerturbationLog(:,:));
    Perturbation_ICU=LBpi_ICU+(UBpi_ICU-LBpi_ICU).*lhs;

    
    %perturbation array Prob Death
    %lb_ProbDeath=log(LBpi_ProbDeath*ones(1,1));
    %ub_ProbDeath=log(UBpi_ProbDeath*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    %PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_ProbDeath-lb_ProbDeath),lb_ProbDeath);
    %Perturbation_ProbDeath(:,:)=exp(PerturbationLog(:,:));
    Perturbation_ProbDeath=LBpi_ProbDeath+(UBpi_ProbDeath-LBpi_ProbDeath).*lhs;
    
    
    %perturbation array presymptomatic period (PresymPeriod)
    lb_Presym=log(LBpi_Presym*ones(1,1));
    ub_Presym=log(UBpi_Presym*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_Presym-lb_Presym),lb_Presym);
    Perturbation_Presym(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_Presym=LBpi_Presym+(UBpi_Presym-LBpi_Presym).*lhs;
    
    
    %perturbation array s0_0
    lb_s0_0=log(LBpi_s0_0*ones(1,1));
    ub_s0_0=log(UBpi_s0_0*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_s0_0 - lb_s0_0),lb_s0_0);
    Perturbation_s0_0(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_s0=LBpi_s0+(UBpi_s0-LBpi_s0).*lhs;
      
    %perturbation array s0_1
    lb_s0_1=log(LBpi_s0_1*ones(1,1));
    ub_s0_1=log(UBpi_s0_1*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_s0_1 - lb_s0_1),lb_s0_1);
    Perturbation_s0_1(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_s0_1=LBpi_s0_1+(UBpi_s0_1-LBpi_s0_1).*lhs;
    
    %perturbation array s0_2
    lb_s0_2=log(LBpi_s0_2*ones(1,1));
    ub_s0_2=log(UBpi_s0_2*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_s0_2 - lb_s0_2),lb_s0_2);
    Perturbation_s0_2(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_s0_2=LBpi_s0_2+(UBpi_s0_2-LBpi_s0_2).*lhs;
    
    %perturbation array s0_3
    lb_s0_3=log(LBpi_s0_3*ones(1,1));
    ub_s0_3=log(UBpi_s0_3*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_s0_3 - lb_s0_3),lb_s0_3);
    Perturbation_s0_3(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_s0_3=LBpi_s0_3+(UBpi_s0_3-LBpi_s0_3).*lhs;
    
    %perturbation array s0_4
    lb_s0_4=log(LBpi_s0_4*ones(1,1));
    ub_s0_4=log(UBpi_s0_4*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_s0_4 - lb_s0_4),lb_s0_4);
    Perturbation_s0_4(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_s0_4=LBpi_s0_4+(UBpi_s0_4-LBpi_s0_4).*lhs;
    
    %perturbation array s0_5
    lb_s0_5=log(LBpi_s0_5*ones(1,1));
    ub_s0_5=log(UBpi_s0_5*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_s0_5 - lb_s0_5),lb_s0_5);
    Perturbation_s0_5(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_s0_5=LBpi_s0_5+(UBpi_s0_5-LBpi_s0_5).*lhs;
    
    %perturbation array s0_6
    lb_s0_6=log(LBpi_s0_6*ones(1,1));
    ub_s0_6=log(UBpi_s0_6*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_s0_6 - lb_s0_6),lb_s0_6);
    Perturbation_s0_6(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_s0_6=LBpi_s0_6+(UBpi_s0_6-LBpi_s0_6).*lhs;

    %perturbation array for fraction FracCrit2
    lb_FracCrit2=log(LBpi_FracCrit2*ones(1,1));
    ub_FracCrit2=log(UBpi_FracCrit2*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_FracCrit2-lb_FracCrit2),lb_FracCrit2);
    Perturbation_FracCrit2(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_FracCrit2=LBpi_FracCrit2+(UBpi_FracCrit2-LBpi_FracCrit2).*lhs;
    
    %perturbation array for fraction FracCrit3
    lb_FracCrit3=log(LBpi_FracCrit3*ones(1,1));
    ub_FracCrit3=log(UBpi_FracCrit3*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_FracCrit3-lb_FracCrit3),lb_FracCrit3);
    Perturbation_FracCrit3(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_FracCrit3=LBpi_FracCrit3+(UBpi_FracCrit3-LBpi_FracCrit3).*lhs;

    %perturbation array for fraction FracCrit4
    lb_FracCrit4=log(LBpi_FracCrit4*ones(1,1));
    ub_FracCrit4=log(UBpi_FracCrit4*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_FracCrit4-lb_FracCrit4),lb_FracCrit4);
    Perturbation_FracCrit4(:,:)=exp(PerturbationLog(:,:));    
    %Perturbation_FracCrit4=LBpi_FracCrit4+(UBpi_FracCrit4-LBpi_FracCrit4).*lhs;

    %perturbation array for fraction FracCrit5
    lb_FracCrit5=log(LBpi_FracCrit5*ones(1,1));
    ub_FracCrit5=log(UBpi_FracCrit5*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_FracCrit5-lb_FracCrit5),lb_FracCrit5);
    Perturbation_FracCrit5(:,:)=exp(PerturbationLog(:,:)); 
    %Perturbation_FracCrit5=LBpi_FracCrit5+(UBpi_FracCrit5-LBpi_FracCrit5).*lhs;

    %perturbation array for transmission parameter be2
    lb_be2=log(LBpi_be2*ones(1,1));
    ub_be2=log(UBpi_be2*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_be2-lb_be2),lb_be2);
    Perturbation_be2(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_be2=LBpi_be2+(UBpi_be2-LBpi_be2).*lhs;    
    
    %perturbation array for transmission parameter b02
    lb_b02=log(LBpi_b02*ones(1,1));
    ub_b02=log(UBpi_b02*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_b02-lb_b02),lb_b02);
    Perturbation_b02(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_b02=LBpi_b02+(UBpi_b02-LBpi_b02).*lhs;    

    
    %perturbation array Prob DeathH
    %lb_ProbDeathH=log(LBpi_ProbDeathH*ones(1,1));
    %ub_ProbDeathH=log(UBpi_ProbDeathH*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    %PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_ProbDeathH-lb_ProbDeathH),lb_ProbDeathH);
    %Perturbation_ProbDeathH(:,:)=exp(PerturbationLog(:,:));
    Perturbation_ProbDeathH=LBpi_ProbDeathH+(UBpi_ProbDeathH-LBpi_ProbDeathH).*lhs;

    %perturbation array n 
    %lb_nHILL=log(LBpi_nHILL*ones(1,1));
    %ub_nHILL=log(UBpi_nHILL*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    %PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_nHILL-lb_nHILL),lb_nHILL);
    %Perturbation_nHILL(:,:)=exp(PerturbationLog(:,:));
    Perturbation_nHILL=LBpi_nHILL+(UBpi_nHILL-LBpi_nHILL).*lhs;
    
    %perturbation array K 
    %lb_KHILL=log(LBpi_KHILL*ones(1,1));
    %ub_KHILL=log(UBpi_KHILL*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    %PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_KHILL-lb_KHILL),lb_KHILL);
    %Perturbation_KHILL(:,:)=exp(PerturbationLog(:,:));
    Perturbation_KHILL=LBpi_KHILL+(UBpi_KHILL-LBpi_KHILL).*lhs;

   
    Perturbation = [Perturbation_be1 Perturbation_b01 Perturbation_b1 Perturbation_b2 Perturbation_b3 Perturbation_FracSev Perturbation_FracCrit Perturbation_FracAsym Perturbation_Incub Perturbation_DurMild Perturbation_DurAsym Perturbation_DurHosp Perturbation_ICU Perturbation_ProbDeath Perturbation_Presym Perturbation_s0_0 Perturbation_s0_1 Perturbation_s0_2 Perturbation_s0_3 Perturbation_s0_4 Perturbation_s0_5 Perturbation_s0_6 Perturbation_FracCrit2 Perturbation_FracCrit3 Perturbation_FracCrit4 Perturbation_FracCrit5 Perturbation_be2 Perturbation_b02 Perturbation_ProbDeathH Perturbation_nHILL Perturbation_KHILL];    
    
    Nominal_I2=zeros(251,Nsample);
    Nominal_I3=zeros(251,Nsample);
    Nominal_D=zeros(251,Nsample);
    
 
    Results{k,1}=Perturbation;
    
    %Perform perturbed simulations and store data
    parfor i=1:Nsample
        piter=pN.*Perturbation(i,:);    
                
        be1 = piter(1);       
        b01 = piter(2);  
        b1 = piter(3);  
        b2 = piter(4);  
        b3 = piter(5);  
        
        be1=be1/10^5;
        b01=b01/10^5;
        b1=b1/10^5;
        b2=b2/10^5;
        b3=b3/10^5;
        
        FracSevere = piter(6);  
        FracCritical = piter(7);  
        FracMild = 1-FracSevere-FracCritical;
        FracAsym = piter(8);  
        
        IncubPeriod = piter(9);  
        DurMildInf = piter(10); 
        DurAsym = piter(11);  
    
        DurHosp = piter(12); 
        TimeICUDeath = piter(13);   
        ProbDeath = piter(14);  
        PresymPeriod = piter(15)*IncubPeriod;           
        CFR=ProbDeath*FracCritical/100; 
            
        a1=1/PresymPeriod; 
        a0=1/(IncubPeriod-PresymPeriod); 

        f=FracAsym;
    
        g0=1/DurAsym;

        g1=(1/DurMildInf)*FracMild;
        p1=(1/DurMildInf)-g1;
    
        p2=(1/DurHosp)*(FracCritical/(FracSevere+FracCritical));
        g2=(1/DurHosp)-p2;
    
        u=(1/TimeICUDeath)*(CFR/FracCritical);
    
        g3=(1/TimeICUDeath)-u;
        
        FracCritical2=piter(23);
        FracCritical3=piter(24);   
        FracCritical4=piter(25); 
        FracCritical5=piter(26);    
        
        be2 = piter(27);       
        b02 = piter(28);      
        
        be2=be2/10^5;
        b02=b02/10^5;
        
        ProbDeathH=piter(29);
        
        n=piter(30);
        K=piter(31);
        
        %how many days after second dose to acquire immunity
        tau_imm = 14; 
        %efficacy first dose
        rho_firstdose = 0.8;  
        %efficacy second dose
        rho_seconddose = 0.95;  
        %days between first and second dose
        tau_first_second_dose = 21; 
        
        pNi = [be1 b01 b1 b2 b3 a0 a1 f g0 g1 p1 g2 p2 g3 u DurHosp FracSevere FracCritical DurMildInf ProbDeath TimeICUDeath FracCritical2 FracCritical3 FracCritical4 FracCritical5 be2 b02 ProbDeathH n K tau_imm rho_firstdose rho_seconddose tau_first_second_dose]; 
        
        s0=[piter(16) piter(17) piter(18) piter(19) piter(20) piter(21) piter(22)];
        s1=1;
        s2=1;
        
        modelsim=@(t,x) ode_covid19_v3_vaccine_umbria(t,x,pNi,Tlock,s0,s1,s2,region_name);
        [T, y]=ode15s(modelsim,time,x0);
        %y(1,10)=NP0;
        %for z=2:length(time)
            %y(z,10)=y(z,4)+y(z,5)+y(z,6)+y(z,7)+y(z,8)+y(z,9)-y(z-1,4)-y(z-1,5)-y(z-1,6)-y(z-1,7)-y(z-1,8)-y(z-1,9);
        %end

             
        msgstr=lastwarn;
        if isequal(msgstr,'')
            Nominal_I2(:,i)=y(:,6);
            Nominal_I3(:,i)=y(:,7);
            Nominal_D(:,i)=y(:,9);  
        else
            Nominal_I2(:,i)=nan(251,1);
            Nominal_I3(:,i)=nan(251,1);
            Nominal_D(:,i)=nan(251,1);
        end
        lastwarn('')
    end
    
    Results{k,2}=Nominal_I2;
    Results{k,3}=Nominal_I3;
    Results{k,4}=Nominal_D;

     
end

delete(poolobj);

toc;

save ('sim_model_covid_Nr=10_Ns=100000_t=250_LINLOG_step1_UMBRIA_I2I3DE_1set_2MAY_OK.mat','-v7.3');
