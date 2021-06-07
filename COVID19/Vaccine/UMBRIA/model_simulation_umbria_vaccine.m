

init_umbria_vaccine
%time=0:1:250;

time_points=1:1:243;
umbria=readtable('data_covid19_umbria_2may.csv');
umbria_matrix=[umbria{:,5},umbria{:,3:4},umbria{:,7}];
  
N = 882000;         %population Umbria
umbria_matrix=(umbria_matrix/N)*10^5;

Nr=10;
for i=1:Nr

        p0=moda_parametri_inizio11step_LINLOG(1,:);
        be1 = p0(1);
        b01=p0(2);
        b1=p0(3);
        b2=p0(4);
        b3=p0(5);
            
        be1=be1/10^5;
        b01=b01/10^5;
        b1=b1/10^5;
        b2=b2/10^5;
        b3=b3/10^5;
        
        FracSevere = p0(6);
        FracCritical = p0(7);
        FracMild = 1-FracSevere-FracCritical;
        FracAsym = p0(8);
        
        IncubPeriod = p0(9);
        DurMildInf = p0(10);
        DurAsym = p0(11);
        
        DurHosp = p0(12);
        TimeICUDeath = p0(13);
        ProbDeath=p0(14);
        PresymPeriod=p0(15)*IncubPeriod; 

        CFR=ProbDeath*FracCritical/100; 
            
        a1=1/PresymPeriod; %presymptomatic period of transmission
        a0=1/(IncubPeriod-PresymPeriod); %true latent period  

        f=FracAsym;
    
        g0=1/DurAsym;
% 
%         g1=(1/DurMildInf)*FracMild;
%         p1=(1/DurMildInf)-g1;
%     
%         p2=(1/DurHosp)*(FracCritical/(FracSevere+FracCritical));
%         g2=(1/DurHosp)-p2;
%     
%         u=(1/TimeICUDeath)*(CFR/FracCritical);
%     
%         g3=(1/TimeICUDeath)-u;
        
        FracCritical2=p0(23);
        FracCritical3=p0(24);
        FracCritical4=p0(25);
        FracCritical5=p0(26);
        be2 = p0(27);       
        b02 = p0(28);   
        be2=be2/10^5;
        b02=b02/10^5;
        
        ProbDeathH=p0(29);
        
        n=p0(30);
        K=p0(31); 
        
        %days after second dose to acquire immunity
        tau_imm = 14; 
        %efficacy of first dose
        rho_firstdose = 0.8;  
        %efficacy of second dose
        rho_seconddose = 0.95;  
        %days between first and second dose
        tau_first_second_dose = 21; 
        
        
        pN = [be1 b01 b1 b2 b3 a0 a1 f g0 DurHosp FracSevere FracCritical DurMildInf ProbDeath TimeICUDeath FracCritical2 FracCritical3 FracCritical4 FracCritical5 be2 b02 ProbDeathH n K tau_imm rho_firstdose rho_seconddose tau_first_second_dose]; 

        s0=[p0(16) p0(17) p0(18) p0(19) p0(20) p0(21) p0(22)];
        s1=1;
        s2=1;
        
        modelsim=@(t,x) ode_covid19_v3_vaccine_umbria(t,x,pN,Tlock,s0,s1,s2,region_name);
        [T, y]=ode15s(modelsim,time,x0);
        
         figure(1)
         for j=1:12
           subplot(4,3,j)  
           plot(time,y(:,j));
           hold on;
           if(j==6)
             scatter(time_points,umbria_matrix(184:end,2),'k');
             hold on
           end
           if(j==7)
             scatter(time_points,umbria_matrix(184:end,3),'k');
             hold on
           end
           if(j==9)
             scatter(time_points,umbria_matrix(184:end,4),'k');
             hold on
           end
         end
 end

figure(2)
plot(time,log(y(:,6)));
hold on;
scatter(time_points,log(umbria_matrix(184:end,2)),'k');

figure(3)
plot(time,log(y(:,7)));
hold on;
scatter(time_points,log(umbria_matrix(184:end,3)),'k');

figure(4)
plot(time,log(y(:,9)));
hold on;
scatter(time_points,log(umbria_matrix(184:end,4)),'k');
 

%R0=10^5*(be/a1 + f*(b0-g0) +(1-f)*1/(p1+g1)*(b1+p1/(p2+g2+u1)*(b2+b3*p2/(u+g3))));
