
 italy_0305=readtable('data_italy_24feb_3may.csv');
 italy_0305_matrix=[italy_0305{:,3:4},italy_0305{:,7}];
  
 N=60e6;
 italy_0305_matrix=(italy_0305_matrix/N)*10^5;