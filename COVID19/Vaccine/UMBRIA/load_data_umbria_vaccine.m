 

umbria=readtable('data_covid19_umbria_2may.csv');
umbria_matrix=[umbria{:,3:4},umbria{:,7}];

N=882000;
umbria_matrix=(umbria_matrix/N)*10^5;