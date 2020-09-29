%ODE model to describe the spread and clinical progression of COVID-19

function dx=ode_covid19_v3_OK(t,x,p,Tlock,s00,s11,s22,region_name)

be = p(1);
b0 = p(2);                
b1 = p(3);
b2 = p(4);
b3 = p(5);
a0 = p(6);
a1 = p(7);
f = p(8);
g0 = p(9);
g1 = p(10);
p1 = p(11);
g2 = p(12);
p2 = p(13);
g3 = p(14);
u = p(15);


[s0,s1,s2]=compute_intervention(t,Tlock,s00,s11,s22,region_name);


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



dS = - ((be * s0)* E1 + (b0*s0)*I0 + (b1*s1)*I1 + (b2*s2)*I2 + (b3*s2)*I3)*S; % + k_e0s * E0;
dE0 = ((be *s0)* E1 + (b0*s0)*I0 + (b1*s1)*I1 + (b2*s2)*I2 + (b3*s2)*I3)*S -a0*E0; % - k_e0s * E0; 
dE1 = a0 * E0 - a1*E1; % -g4*E1;
dI0 = f * a1 * E1 - g0*I0;
dI1 = (1-f) * a1 * E1 - g1*I1 - p1*I1;
dI2 = p1 * I1 - g2*I2 - p2*I2;  % -u2*I2;
dI3 = p2 * I2 - g3 *I3 - u * I3 ;
dR = g0*I0 + g1*I1 + g2*I2 + g3*I3; % +g4 * E1;
dD = u * I3; % + u2 * I2;




dx=[dS;dE0;dE1;dI0;dI1;dI2;dI3;dR;dD];
