observablesName={'I2','I3','DE'};

Results_new=error_tot_LINLOG_10Nr_1step_ITALY';

%ITALIA - threshold list

lowThr_1step=[1750,266,1406]; 
%lowThr_2step=[1720,255,1390]; 
%lowThr_3step=[1650,230,1280]; 
%lowThr_4step=[1430,200,1130]; 
%lowThr_5step=[1250,100,1050]; 
%lowThr_6step=[900,70,850]; 
%lowThr_7step=[680,55,700]; 
%lowThr_8step=[600,55,300]; 



%UMBRIA - threshold list

%lowThr_1step=[572,180,235];
%lowThr_2step=[565,172,229];
%lowThr_3step=[525,160,219];
%lowThr_4step=[320,115,160];
%lowThr_5step=[200,91,91];
%lowThr_6step=[176,63,71];
%lowThr_7step=[165,48,57];  
%lowThr_8step=[150,41,44];


highThr=6000000*ones(1,3);

indexes=upperLowerSet_error_Nr(Results_new,Nr,Nsample,observablesName,lowThr_1step,highThr);

all_intersection={};
 for k=1:Nr
   k
  intersection=combinazioni(indexes(2*k-1,:));      
  all_intersection=[all_intersection;intersection];
 end

 
 
MIRIT=[];
CloudCondUUUT=[];
CloudCondLLLT=[];


xbin_p_1=cell(1,10);
ks_p_1=cell(1,10);

for k=1:Nr
    k
    NCloudUUU=0; 
    NCloudLLL=0;
    
   
    for i=1:length(all_intersection(k,:))
        if strcmp(all_intersection{2*k,i},'TTT')
            NCloudUUU=length(all_intersection{2*k-1,i});
            index_1=i; 
            break;
        end
    end
       
    for i=1:length(all_intersection(k,:))
        if strcmp(all_intersection{2*k,i},'FFF')
            NCloudLLL=length(all_intersection{2*k-1,i});
            index_2=i;
            break;
        end
    end   
 
    Ncloud=NCloudLLL;
    
   
    %pN=mode_parameters_inizio2step_LINLOG(k,:);
    param_intersection=pN.*Results{k,1}(all_intersection{2*k-1,index_2},:);
 
    [ks_LLL, xbin_LLL]=mvksdensity(param_intersection,param_intersection,'Bandwidth',0.8);    
    xbin_p_1{1,k}=xbin_LLL;
    ks_p_1{1,k}=ks_LLL;   
    
            
    [ks_LLL_m, ks_LLL_I]=max(ks_LLL);
    CloudCondLLL=xbin_LLL(ks_LLL_I,:);
    CloudCondLLLT=[CloudCondLLLT;CloudCondLLL];

end

mode_parameters_inizio2step_LINLOG = CloudCondLLLT;
save('mode_parameter_2step_CRC_italy.mat','mode_parameters_inizio2step_LINLOG');
