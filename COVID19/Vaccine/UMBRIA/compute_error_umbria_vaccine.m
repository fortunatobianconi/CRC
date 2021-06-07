
 load_data_umbria_vaccine
 
 Nr=10;
 Nsample=100000;

error_tot_LINLOG_10Nr_1step_UMBRIA=cell(3,Nr);
 for k=1:Nr
     k
     for i=1:3
         for j=1:Nsample
             error_tot_LINLOG_10Nr_1step_UMBRIA{i,k}(1,j)= (sum(abs(umbria_matrix(184:416,i) - Results{k,i+1}(1:233,j)))); 
         end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
     end                                           
 end
 
 
%normalized error 
error_tot_LINLOG_10Nr_1step_NORM_UMBRIA=cell(3,Nr);
  for k=1:Nr
      k
      for i=1:3
          for j=1:Nsample
              error_tot_LINLOG_10Nr_1step_NORM_UMBRIA{i,k}(1,j)= 1/232* sum(abs(umbria_matrix(184:416,i) - Results{k,i+1}(1:233,j))./max(1,umbria_matrix(184:416,i)));
          end
      end
  end

 save('error_1step_LINLOG_10Nr_100000_UMBRIA_I2I3DE_1SET_2may.mat','error_tot_LINLOG_10Nr_1step_UMBRIA');
 save('error_1step_LINLOG_10Nr_100000_NORM_UMBRIA_I2I3DE_1SET_2may.mat','error_tot_LINLOG_10Nr_1step_NORM_UMBRIA');
  
 