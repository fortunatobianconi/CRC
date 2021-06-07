
observablesName={'I2','I3','DE'};

Results_new=error_tot_LINLOG_10Nr_1step_UMBRIA';

%1 september 2 may
lowThr_1step=[4300,720,7300];
%lowThr_2step=[3700,690,6800];
%lowThr_3step=[3500,600,6500];
%lowThr_4step=[3300,580,6400];
%lowThr_5step=[3000,570,6300];
%lowThr_6step=[2300,430,5000];
%lowThr_7step=[1800,350,4400];
%lowThr_8step=[1300,250,2000];
%lowThr_9step=[1200,200,1800];
%lowThr_10step=[800,170,1000];

highThr=6000000*ones(1,4);

indexes=upperLowerSet_error_Nr(Results_new,Nr,Nsample,observablesName,lowThr_1step,highThr);

all_intersection={};
 for k=1:Nr
   
  intersection=combinazioni(indexes(2*k-1,:));      
  all_intersection=[all_intersection;intersection];
 end

 
 
MIRIT=[];
CloudCondUUUT=[];
CloudCondLLLT=[];


xbin_p_1=cell(1,1);
ks_p_1=cell(1,1);

for k=1:Nr
   
    NCloudUUU=0; 
    NCloudLLL=0;
    
   
    for i=1:length(all_intersection(k,:))
        if strcmp(all_intersection{2*k,i},'TTT') %high
            NCloudUUU=length(all_intersection{2*k-1,i});
            index_1=i; 
            break;
        end
    end
       
    for i=1:length(all_intersection(k,:))
        if strcmp(all_intersection{2*k,i},'FFF')  %low
            NCloudLLL=length(all_intersection{2*k-1,i});
            index_2=i;
            break;
        end
    end   
 


    Ncloud=NCloudLLL;  
    
    %Uncomment from second step onwards
    %pN=moda_parametri_inizio2step_LINLOG(k,:);
    param_intersection=pN.*Results{k,1}(all_intersection{2*k-1,index_2},:);
 
    [ks_LLL, xbin_LLL]=mvksdensity(param_intersection,param_intersection,'Bandwidth',0.8);    
    xbin_p_1{1,k}=xbin_LLL;
    ks_p_1{1,k}=ks_LLL;

    %plot univariate pdf for each parameter
%     for ip=1:size(pN,2)
%   
%         parameter_axis=min(pN(ip).*Results{k,1}(:,ip)):(max(pN(ip).*Results{k,1}(:,ip)) - min(pN(ip).*Results{k,1}(:,ip)))/Ncloud:max(pN(ip).*(Results{k,1}(:,ip)));
%         param_univ = pN(ip).*Results{k,1}(all_intersection{2*k-1,index_2},ip);
%         [ks_LLL_univ xbin_LLL_univ]=ksdensity(param_univ,parameter_axis);        
%         
%         figure
%         plot(xbin_LLL_univ,ks_LLL_univ,'b');
%     end
     
            
    [ks_LLL_m, ks_LLL_I]=max(ks_LLL);
    CloudCondLLL=xbin_LLL(ks_LLL_I,:);
    CloudCondLLLT=[CloudCondLLLT;CloudCondLLL];

end


mode_parameters_2step_LINLOG = CloudCondLLLT;

save('mode_parameters_2step_LINLOG_Nr10_Ns100000_UMBRIA_I2I3DE_multivariate_1SET_2may.mat','mode_parameters_2step_LINLOG');
save('intersection_1step_LINLOG__Nr10_I2I3DE_UMBRIA_1SET_2may.mat','all_intersection');
save('pdf_param_1step_LINLOG__Nr10_I2I3DE_UMBRIA_1SET_2may.mat','ks_p_1','xbin_p_1');