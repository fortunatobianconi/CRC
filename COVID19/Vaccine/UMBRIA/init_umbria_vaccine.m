region_name='Umbria_second_third_wave';

%define time
time = 0:1:250;  

%define initial conditions
N = 882000;         %population Umbria

%1 september
E00 = (400/N)*10^5;
S0 = 10^5-E00;

%1 september 
x0 = [S0 E00 (300/N)*10^5 (144/N)*10^5 (144/N)*10^5 (8/N)*10^5 (2/N)*10^5 (1452/N)*10^5 (80/N)*10^5 0 0 0]; 

%14 September: 14days school opening
%19 October:49days DPCM in Umbria
%11 November: 72days orange zone
%6 December: 97days yellow zone
%24 December: 115days red zone 
%7 January: 129days end of Christmas holidays 
%8 February: 160days red zone in the Province of Perugia
%22 March: 202days yellow zone and school reopening
%26 April: 237days yellow zone
%22 May: 265days reopening of gym, shopping centers, ...
%21 June: 293days no curfew

Tlock=[14 49 72 97 129 160 202 ]; 
