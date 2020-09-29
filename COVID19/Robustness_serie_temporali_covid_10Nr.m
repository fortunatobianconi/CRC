%SCRIPT FOR RUNNING CRC USING THE ITALIAN DATA 


tic;

init_italy
% Load parameters (set to 1 or mode previous step)
%load('mode_parameter_2step_CRC_italy.mat');
pN=ones(1,21);

% Generate poerturbation matrix usign latin hypercube sampling
Nsample=100000; % number of parameter samples
Nr=10; % number of realizations

Results=cell(Nr,4);

% Performs perturbed simulations and store data
poolobj=parpool(4);
for k=1:Nr
    k
    
    %pN=moda_parametri_inizio2step_LINLOG(k,:);
    
    % Latin hypercube perturbation for transmission rate (be,b0,b1,b2,b3)    
    LBpi_be=0.001/pN(1); % Lower bound   
    UBpi_be=3/pN(1);   % Upper bound  
 
    LBpi_b0=0.001/pN(2); % Lower bound   
    UBpi_b0=3/pN(2);   % Upper bound  
 
    LBpi_b1=0.001/pN(3); % Lower bound    
    UBpi_b1=3/pN(3);   % Upper bound  
 
    LBpi_b2=0.001/pN(4); % Lower bound    
    UBpi_b2=3/pN(4);   % Upper bound  
 
    LBpi_b3=0.001/pN(5); % Lower bound  
    UBpi_b3=3/pN(5);   % Upper bound  
 
    % bounds for fractions (FracSevere, FracCritical, FracMild)
    LBpi_FracSev=0.01/pN(6);  
    UBpi_FracSev=0.4/pN(6);   
 
    LBpi_FracCrit=0.01/pN(7);   
    UBpi_FracCrit=0.3/pN(7);    
 
    %bounds for fractions FracAsym
    LBpi_FracAsym = 0.1/pN(8);   
    UBpi_FracAsym = 0.6/pN(8);   
 
    %bounds for incubation period, mild infection period and asymptomatic
    %infection period (IncubPeriod, DurMildInf,DurAsym)
    LBpi_Incub=4/pN(9);    
    UBpi_Incub=14/pN(9);   
 
    LBpi_DurMild=2/pN(10);  
    UBpi_DurMild=80/pN(10);  
 
    LBpi_DurAsym=2/pN(11);   
    UBpi_DurAsym=30/pN(11);    
 
    %bounds for duration severe infection (DurHosp)
    LBpi_DurHosp=2/pN(12);   
    UBpi_DurHosp=90/pN(12);   
 
    %bounds for ICU stay (TimeICUDeath)
    LBpi_ICU=2/pN(13);   
    UBpi_ICU=70/pN(13);   
 
    %bounds for ProbDeath
    LBpi_ProbDeath=1/pN(14);   
    UBpi_ProbDeath=90/pN(14);  
 
    %bounds for PreSymPeriod
    LBpi_Presym=0.5/pN(15);  
    UBpi_Presym=0.9/pN(15);     
    
    %bounds for s0_0
    LBpi_s0_0=0.5/pN(16);  
    UBpi_s0_0=0.9/pN(16);  
    
    %bounds for s0_1
    LBpi_s0_1=0.4/pN(17);  
    UBpi_s0_1=0.9/pN(17);  
    
    %bounds for s0_2
    LBpi_s0_2=0.3/pN(18);  
    UBpi_s0_2=0.7/pN(18);  
    
    %bounds for s0_3
    LBpi_s0_3=0.05/pN(19);  
    UBpi_s0_3=0.5/pN(19);
 
    %bounds for s1
    LBpi_s1=0.1/pN(20);  
    UBpi_s1=0.9/pN(20);  
    
    %bounds for s2
    LBpi_s2=0.1/pN(21);  
    UBpi_s2=0.9/pN(21);  


    %perturbation array for transmission rate parameters (be,b0,b1,b2,b3)
    lb_be=log(LBpi_be*ones(1,1));
    ub_be=log(UBpi_be*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_be-lb_be),lb_be);
    Perturbation_be(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_be=LBpi_be+(UBpi_be-LBpi_be).*lhs;
    
    
    lb_b0=log(LBpi_b0*ones(1,1));
    ub_b0=log(UBpi_b0*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_b0-lb_b0),lb_b0);
    Perturbation_b0(:,:)=exp(PerturbationLog(:,:));
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
    
    
    %perturbation array for fraction parameters (FracSevere, FracCritical, FracMild, FracAsym)
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
    %Perturbation_s0=LBpi_s0+(UBpi_s0-LBpi_s0).*lhs;
    
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
     
    %perturbation array s1_0
    lb_s1=log(LBpi_s1*ones(1,1));
    ub_s1=log(UBpi_s1*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_s1-lb_s1),lb_s1);
    Perturbation_s1(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_s1=LBpi_s1+(UBpi_s1-LBpi_s1).*lhs;
    
   
    %perturbation array s2
    lb_s2=log(LBpi_s2*ones(1,1));
    ub_s2=log(UBpi_s2*ones(1,1));
    lhs=lhsdesign(Nsample,1);
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub_s2-lb_s2),lb_s2);
    Perturbation_s2(:,:)=exp(PerturbationLog(:,:));
    %Perturbation_s2=LBpi_s2+(UBpi_s2-LBpi_s2).*lhs;

    
    Perturbation = [Perturbation_be Perturbation_b0 Perturbation_b1 Perturbation_b2 Perturbation_b3 Perturbation_FracSev Perturbation_FracCrit Perturbation_FracAsym Perturbation_Incub Perturbation_DurMild Perturbation_DurAsym Perturbation_DurHosp Perturbation_ICU Perturbation_ProbDeath Perturbation_Presym Perturbation_s0_0 Perturbation_s0_1 Perturbation_s0_2 Perturbation_s0_3 Perturbation_s1 Perturbation_s2];
    
    Nominal_I2=zeros(111,Nsample);
    Nominal_I3=zeros(111,Nsample);
    Nominal_D=zeros(111,Nsample);
    
 
    Results{k,1}=Perturbation;

    parfor i=1:Nsample
        piter=pN.*Perturbation(i,:);    
        
        be = piter(1);
        b0=piter(2);
        b1=piter(3);
        b2=piter(4);
        b3=piter(5);
        
        be=be/10^5;
        b0=b0/10^5;
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
        ProbDeath=piter(14);
        PresymPeriod=piter(15)*IncubPeriod;
        
        CFR=ProbDeath*FracCritical/100; 
            
        a1=1/PresymPeriod; %presymptomatic period of transmission
        a0=1/(IncubPeriod-PresymPeriod); %true latent period  

        f=FracAsym;
    
        g0=1/DurAsym;

        g1=(1/DurMildInf)*FracMild;
        p1=(1/DurMildInf)-g1;
    
        p2=(1/DurHosp)*(FracCritical/(FracSevere+FracCritical));
        g2=(1/DurHosp)-p2;
    
        u=(1/TimeICUDeath)*(CFR/FracCritical);
    
        g3=(1/TimeICUDeath)-u;

        pNi = [be b0 b1 b2 b3 a0 a1 f g0 g1 p1 g2 p2 g3 u]; 
        
        s0 = [piter(16) piter(17) piter(18) piter(19)];  
        s1=piter(20);  
        s2=piter(21); 
        
        modelsim=@(t,x) ode_covid19_v3_final_perCRC(t,x,pNi,Tlock,s0,s1,s2,region_name);
        [T, y]=ode15s(modelsim,time,x0);
             
        msgstr=lastwarn;
        if isequal(msgstr,'')           
            %Nominal_I1(:,i)=y(:,5);
            Nominal_I2(:,i)=y(:,6);
            Nominal_I3(:,i)=y(:,7);
            %Nominal_R(:,i)=y(:,8);
            Nominal_D(:,i)=y(:,9);    
            
            
        else
            
            %Nominal_I1(:,i)=nan(111,1);
            Nominal_I2(:,i)=nan(111,1);
            Nominal_I3(:,i)=nan(111,1);
            %Nominal_R(:,i)=nan(111,1);
            Nominal_D(:,i)=nan(111,1);
                        
            
        end
        lastwarn('')
    end
    %Results{k,2}=Nominal_I1;
    Results{k,2}=Nominal_I2;
    Results{k,3}=Nominal_I3;
    %Results{k,5}=Nominal_R;
    Results{k,4}=Nominal_D;

     
end

delete(poolobj);

toc;

save ('results_italy_step1.mat','-v7.3');