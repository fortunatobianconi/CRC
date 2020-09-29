region_name='Umbria';

%define time
time = 0:1:80;

%define initial conditions
N = 882000;         %population Umbria

E00 = (1/N)*10^5;
S0 = 10^5-E00;
x0 = [S0 E00 0 0 0 0 0 0 0];
NP0 = 0;

        
%start of the epidemic: 21 February
%after 3 days (24 February) PPE
%after 13 days (5 March): school closure
%after 16 days (8 March): social distancing + no events
%after 18 days (10 March): total lock-down


Tlock=[3 13 18];

