
region_name='Italy';

%define time
time = 0:1:110;

%define initial conditions
N = 60e6;         %population Italy

E00 = (1/N)*10^5;
S0 = 10^5-E00;
x0 = [S0 E00 0 0 0 0 0 0 0];
NP0=0;


%epidemics start: 27 January
%after 25 days (21 February): red areas in Lombardy and Veneto
%aftre 28 days (24 February): school closures in Emilia-Romagna, Friuli-Venezia Giulia,
%Liguria, Lombardia, Marche, Piemonte, Veneto
%after 25 days (21 February) PPE
%after 38 days (5 March): school closure in the rest of Italy
%after 41 days (8 March): lockdown Norther of Italy + social distancing + banning of public events 
%after 43 days (10 March): total lock-down

Tlock=[25 28 38 43];
    
