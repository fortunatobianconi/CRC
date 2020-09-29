 

umbria_0305=readtable('data_covid19_umbria_3may.csv');
umbria_0305_matrix=[umbria_0305{:,3:4},umbria_0305{:,7}];

N=882000; %Umbria population
umbria_0305_matrix=(umbria_0305_matrix/N)*10^5;