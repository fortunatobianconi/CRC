%  For each output variable in input, it returns the set of samples of the parameter vector for which
%  the corresponding distance function has a value under the defined tolerance. It takes in input:
%  - Distance: the measured distance functions
%  - Nr: number of realizations
%  - NSample: number of samples of the parameter vector
%  - proteinName: name of measured proteins
%  - lowThr or highThr: specify the maximum or minimum tolerance for each distance function. It is possible to
%    choose a different threshold for each variable


function results=upperLowerSet_Nr(Distance,Nr,Nsample,proteinName,lowThr,highThr)

% final results
results=cell(Nr,2*(size(Distance,2)));

% remove possible negative values
for k=1:Nr  
 for j=1:size(Distance,2)
  for i=1:length(Distance{k,1})
     if(Distance{k,j}(i)<0) 
         Distance{k,j}(i)=NaN;
     end
  end
 end
end


% for each realization
for k=1:Nr
   for s=1:(size(Distance,2))
       high = genvarname(strcat('t_',proteinName{1,s},'_Upper'));  % High values of distance functions
       Upper=[];
       
       low = genvarname(strcat('t_',proteinName{1,s},'_Lower'));  % Low values of distance functions
       Lower=[];
       
      % for each observable, scanning of the results 
      for t=1:(length(Distance{k,s}))
          
          if isempty(highThr)
             if (Distance{k,s}(t)<lowThr(s))
                Lower=[Lower t];
             end
          elseif isempty(lowThr)
               if (Distance{k,s}(t)>highThr(s))
                Upper=[Upper t];
               end    
          else 
              if (Distance{k,s}(t)<lowThr(s))
                Lower=[Lower t];
              end
              if (Distance{k,s}(t)>highThr(s))
                Upper=[Upper t];
               end  
          end 
      end
      
      if isempty(Upper)
          Upper=[0];
      end
      if isempty(Lower)
          Lower=[0];
      end
      
     % results are bounded in a single cell array 
     eval([high ' = [Upper];']);
     eval([low ' = [Lower];']);
     results{2*k-1,2*s-1}=Upper;    
     results{2*k,2*s-1}='Upper';
     results{2*k-1,2*s}=Lower;         
     results{2*k,2*s}='Lower';
   end
end
    
    
