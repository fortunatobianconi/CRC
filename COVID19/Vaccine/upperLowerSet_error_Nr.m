function index=upperLowerSet_Nr(error,Nr,Nsample,proteinName,lowThr,highThr)

index=cell(2*Nr,2*(size(error,2)));

for k=1:Nr  
 for j=1:size(error,2)
  for i=1:length(error{k,1})
     if(error{k,j}(i)<0) %|| Results{k,j}(i)>thresholds(j-1))
       error{k,j}(i)=NaN;
     end
  end
 end
end

% Results_sorted={};
% Results_sorted_indexes={};
% for i=1:length(Results)
%     [Results_sorted{1,i} Results_sorted_indexes{1,i}]=sort(Results{1,i});
% end

for k=1:Nr
    k
  
   for s=1:(size(error,2))
       v = genvarname(strcat('t_',proteinName{1,s},'_Upper'));  %True
       Upper=[];
       
       w = genvarname(strcat('t_',proteinName{1,s},'_Lower'));  %False
       Lower=[];
       
       
      for t=1:(length(error{k,s}))
         if (error{k,s}(t)<lowThr(s)) %&& Results{k,s}(t)>(lowThr))
            Lower=[Lower t];
         else  if (error{k,s}(t)>highThr(s))%lowThr*mean(Results{k,s}(~isnan(Results{k,s})
            Upper=[Upper t];
         %end else if
             end %end if
         end  %end else if
      end %end for interno
      if isempty(Upper)
          Upper=0;
      end
       if isempty(Lower)
          Lower=0;
      end
     eval([v ' = [Upper];']);
     eval([w ' = [Lower];']);
     index{2*k-1,2*s-1}=Upper;        %UPPER             
     index{2*k,2*s-1}='Upper';
     index{2*k-1,2*s}=Lower;          %LOWER
     index{2*k,2*s}='Lower';
   end
end
    
    
