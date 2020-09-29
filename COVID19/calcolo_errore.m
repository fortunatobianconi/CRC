
 load_data_italy
 Nr=10;
 Nsample=100000;
 

error_tot_LINLOG_10Nr_1step_ITALY=cell(3,Nr);
 for k=1:Nr
     k
     for i=1:3
         for j=1:Nsample
             error_tot_LINLOG_10Nr_1step_ITALY{i,k}(1,j)= (sum(abs(italy_0305_matrix(:,i) - Results{k,i+1}(28:97,j))));                                                                                                                      
         end
     end
 end
 
%NORMALIZED ERROR
error_tot_LINLOG_10Nr_1step_NORM_ITALY=cell(3,Nr);
  for k=1:Nr
      k
      for i=1:3
          for j=1:Nsample
              error_tot_LINLOG_10Nr_1step_NORM_ITALY{i,k}(1,j)= 1/69 * sum(abs(italy_0305_matrix(:,i) - Results{k,i+1}(28:97,j))./max(1,italy_0305_matrix(:,i)));
          end
      end
  end
  
 save('error_1step_LINLOG_10Nr_100000_ITALY_I2I3DE_3may.mat','error_tot_LINLOG_10Nr_1step_ITALY');
 save('error_1step_LINLOG_10Nr_100000_NORM_ITALY_I2I3DE_3may.mat','error_tot_LINLOG_10Nr_1step_NORM_ITALY');
  
  