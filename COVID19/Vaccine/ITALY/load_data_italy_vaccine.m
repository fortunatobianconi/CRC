
italy=readtable('data_italy_24feb_2may.csv');
italy_matrix=[italy{:,3:4},italy{:,7}];

N=60*10^6;
italy_matrix=(italy_matrix/N)*10^5;