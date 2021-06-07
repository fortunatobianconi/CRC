region_name='Italy_second_third_wave';

%define time
time = 0:1:250;  

%define initial conditions
N = 60*10^6;         %population Italia


%1 September
x0 = [10^5-(30000/N)*10^5 (30000/N)*10^5 (20000/N)*10^5 (15000/N)*10^5 (26271/N)*10^5 (1437/N)*10^5 (109/N)*10^5 (208201/N)*10^5 (35497/N)*10^5 0 0 0]; 

%NP0 = (1326/N)*10^5;

%14 September: 14days school reopening
%6 November: 68days color coded system
%24 December: 116days red zone in the entire country
%7 January: 130days  end of Christmas holidays
%15 March: 197days removal of yellow zone
%26 April: 237days return to yellow zone
%24 May: 265days gym opening, shopping centers, ...
%21 June: 293days no curfew
Tlock=[14 68 116 130 197 237]; 

