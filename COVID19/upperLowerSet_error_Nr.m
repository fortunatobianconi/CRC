function index=upperLowerSet_Nr(error,Nr,Nsample,proteinName,lowThr,highThr)

%load('Binedges_err_aree.mat');
%load('ks_y_err_aree.mat');


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
%    for i=1:length(proteinName)
%        v = genvarname(strcat('t_',proteinName{1,i},'_Upper'));
%        eval([v ' = []'])
%        
%        w = genvarname(strcat('t_',proteinName{1,i},'_Lower'));
%        eval([w ' = []'])    
%    end  
%    
   for s=1:(size(error,2))
       v = genvarname(strcat('t_',proteinName{1,s},'_Upper'));  %True
       Upper=[];
       
       w = genvarname(strcat('t_',proteinName{1,s},'_Lower'));  %False
       Lower=[];
       
       %BinEdges=min(Results{k,s}):(max(Results{k,s})-min(Results{k,s}))/(Nsample):max(Results{k,s});
       %[ks_y xbin_max]=ksdensity(Results{k,s},BinEdges);
       
           
       %percentileLow=prctile(Results{1,s},lowThr);
       %percentileHigh=prctile(Results{1,s},highThr);
       %percentileLow=lowThr*mean(Results{1,s}(~isnan(Results{1,s})));
       %percentileHigh=highThr*mean(Results{1,s}(~isnan(Results{1,s})));
       %percentileLow=mean(Results{1,s}(~isnan(Results{1,s})))-(std(Results{1,s}(~isnan(Results{1,s})))/mean(Results{1,s}(~isnan(Results{1,s}))));
       %percentileHigh=mean(Results{1,s}(~isnan(Results{1,s})))+(std(Results{1,s}(~isnan(Results{1,s})))/mean(Results{1,s}(~isnan(Results{1,s}))));
       %percentileLow=centralValue(s)-lowThr(s);
       %percentileHigh=centralValue(s)+highThr(s);
       
      for t=1:(length(error{k,s}))
        if (error{k,s}(t)<lowThr(s)) %&& Results{k,s}(t)>(lowThr))
        %if (t<=lowThr) 
        %BinEdges_AMPK_Upper=[BinEdges_AMPK_Upper Results{k,2}(t)];
         %a=[a Results_sorted_indexes{1,s}(t)]; % t è l'indice corrispondente ad un valore alto dell'area
          Lower=[Lower t];
        else  if (error{k,s}(t)>highThr(s))%lowThr*mean(Results{k,s}(~isnan(Results{k,s})
         %BinEdges_AMPK_Lower=[BinEdges_AMPK_Lower Results{k,2}(t)];
         %b=[b Results_sorted_indexes{1,s}(t)];
          Upper=[Upper t];
         %end else if
        end %end if
      end  %end else if
      end %end for interno
      if isempty(Upper)
          Upper=[0];
      end
       if isempty(Lower)
          Lower=[0];
      end
     eval([v ' = [Upper];']);
     eval([w ' = [Lower];']);
     index{2*k-1,2*s-1}=Upper;        %UPPER             %index=[index a b];
     index{2*k,2*s-1}='Upper';
     index{2*k-1,2*s}=Lower;          %LOWER
     index{2*k,2*s}='Lower';
   end
end
    
    
