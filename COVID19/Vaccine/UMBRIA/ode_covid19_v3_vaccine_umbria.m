%ODE model to describe the spread and clinical progression of COVID-19 in
%Umbria

function dx=ode_covid19_v3_vaccine_umbria(t,x,p,Tlock,s00,s11,s22,region_name)

be1 = p(1);
b01 = p(2);     
b1 = p(3);
b2 = p(4);
b3 = p(5);
a0 = p(6);
a1 = p(7);
f = p(8);
g0 = p(9);
% g1 = p(10);
% p1 = p(11);
% g2 = p(12);
% p2 = p(13);
% g3 = p(14);
% u = p(15);
DurHosp = p(10);
FracSevere = p(11);
FracCritical = p(12);
DurMildInf=p(13);
ProbDeath=p(14);
TimeICUDeath=p(15);
FracCritical2=p(16);
FracCritical3=p(17);
FracCritical4=p(18);
FracCritical5=p(19);

be2=p(20);
b02=p(21);  
ProbDeathH=p(22);

n=p(23);
K=p(24);

tau_imm=p(25);
rho_firstdose=p(26);
rho_seconddose=p(27);
tau_first_second_dose=p(28);


%110 Umbria
be=be1*(t<=110) + be2*(t>110);
b0=b01*(t<=110) + b02*(t>110);

Critical=FracCritical*(t<=35) + FracCritical2*(t>35)*(t<=83) + FracCritical3*(t>83)*(t<=152) + FracCritical4* (t>152)*(t<200)+FracCritical5* (t>=200) ;
 
u1=(1/DurHosp)*((ProbDeathH*FracSevere)/100)/FracSevere;
p2 = (1/DurHosp)*(Critical/(FracSevere+Critical));
g2=(1/DurHosp)-p2 -u1;
g1=(1/DurMildInf)*(1-FracSevere-Critical);
p1=(1/DurMildInf)-g1;
CFR=(ProbDeath*Critical/100); 
u=(1/TimeICUDeath)*(CFR/Critical);
g3=(1/TimeICUDeath)-u; 
 

[s0,s1,s2]=compute_intervention_vaccine(t,Tlock,s00,s11,s22,region_name);

%susceptibles individuals (not infected)
S = x(1);
%exposed individuals, infected but no symptoms and no transmission
E0 = x(2);
%exposed individuals, infected but no symptoms and can transmit
E1 = x(3);
%infected individuals but asymptomatic infection
I0 = x(4);
%infected individuals, mild infection
I1 = x(5);
%infected individuals, severe infection
I2 = x(6);
%infected individuals, critical infection
I3 = x(7);
%recovered
R = x(8);
%dead
D = x(9);
%vaccinated first dose
V1 = x(10);
%vaccinated second dose
V2 = x(11);
%immunized
Im = x(12);

%vaccination rate
%8k: 0.0091
%12.5k: 0.0142
%18k: 0.0204
eta = 0*(t<118) + (461/882000) * (t>=118) * (t<154) + (703/882000) * (t>=154) * (t<181) + (2390/882000) * (t>=181)* (t<212) + (3617/882000)* (t>=212)*(t<257) + (18000/882000) * (t>=257); %*(t<500) + (25000/882000)*(t>=500);

%duration of waning immunity
delta=1/240;

Inf = (be * s0)* E1 + (b0*s0)*I0 + (b1*s1)*I1 + (b2*s2)*I2 + (b3*s2)*I3;
dS = - (Inf)*S - eta*S +delta*R + delta*Im  ;
dE0 = (Inf)*S -a0*E0 + (1-rho_firstdose)*Inf*V1 + (1-rho_seconddose)*Inf*V2;
dE1 = a0 * E0 - a1*E1; 
dI0 = f * a1 * E1 - g0*I0;
dI1 = (1-f) * a1 * E1 - g1*I1 - p1*I1;
dI2 = p1 * I1 - g2*I2 - p2*I2 -u1*I2 + 1/(1+(Im/K)^n)*(t>118);
dI3 = p2 * I2 - g3 *I3 - u * I3 ;
dR = g0*I0 + g1*I1 + g2*I2 + g3*I3 -delta*R;
dD = u * I3 + u1 * I2;
dV1 = eta*S - V1/tau_first_second_dose - (1-rho_firstdose)*(Inf)*(V1);
dV2 = V1/tau_first_second_dose - V2/tau_imm - (1-rho_seconddose)*(Inf)*(V2);
dIm = V2/tau_imm -delta*Im;

dx=[dS;dE0;dE1;dI0;dI1;dI2;dI3;dR;dD;dV1;dV2;dIm];