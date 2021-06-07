observablesName={'I2','I3','DE'};

Results_new=error_tot_LINLOG_10Nr_1step_ITALY';

%1 september 2 may
lowThr_1step=[4700,550,7600];
%lowThr_2step=[4500,490,7100];
%lowThr_3step=[4100,460,6700];
%lowThr_4step=[3500,420,6200];
%lowThr_5step=[3000,370,5500];
%lowThr_6step=[2400,290,4000];
%lowThr_7step=[2000,250,3600];
%lowThr_8step=[1600,200,3100];
%lowThr_9step=[1500,150,2000];
%lowThr_10step=[1100,130,1500];
%lowThr_11step=[1100,100,1100];

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
        if strcmp(all_intersection{2*k,i},'FFF') %low
            NCloudLLL=length(all_intersection{2*k-1,i});
            index_2=i;
            break;
        end
    end   
 

    Ncloud=NCloudLLL;
    
    %Uncomment from second step onwards
    %pN=mode_parameters_2step_LINLOG(k,:);
    param_intersection=pN.*Results{k,1}(all_intersection{2*k-1,index_2},:);
        
    [ks_LLL, xbin_LLL]=mvksdensity(param_intersection,param_intersection,'Bandwidth',0.8);    
    xbin_p_1{1,k}=xbin_LLL;
    ks_p_1{1,k}=ks_LLL;    
    
            
    [ks_LLL_m, ks_LLL_I]=max(ks_LLL);
    CloudCondLLL=xbin_LLL(ks_LLL_I,:);
    CloudCondLLLT=[CloudCondLLLT;CloudCondLLL];

end


mode_parameters_2step_LINLOG = CloudCondLLLT;

save('mode_parameters_2step_LINLOG_Nr10_Ns100000_ITALY_I2I3DE_multivariate_1SET_2may.mat','mode_parameters_2step_LINLOG');
save('intersection_2step_LINLOG__Nr10_I2I3DE_ITALY_1SET_2may.mat','all_intersection');
save('pdf_param_2step_LINLOG__Nr10_I2I3DE_ITALY_1SET_2may.mat','ks_p_1','xbin_p_1');