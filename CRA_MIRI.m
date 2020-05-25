% Conditional robustness analysis: script for the computation of MIRI 
%(Moment Independent Robustness Indicator) for each parameter

clear all
close all

% load final mode parameters obtained in the calibration toward low values of the
% distance functions
load('parameters_low.mat')
p_low=parameters_low;

% load final mode parameters obtained in the calibration toward high values of the
% distance functions 
load('parameters_high.mat')
p_high=parameters_high;

% load measures
load('MeasuresLOW.mat');
load('MeasuresHIGH.mat');

% join perturbations of both calibration procedures
for k=1:Nr
Results{k,1}=[PerturbationLOW{k,1};PerturbationHIGH{k,1}];
end

% load the region obtained at the end of each calibration 
load('all_region_high.mat');
load('all_region_low.mat');

MIRIT=[];


% parameter density estimation conditioned to both regions
for k=1:Nr
    
    NCloudU=0; 
    NCloudL=0;
    
    for i=1:length(all_intersection_high(k,:))
        if strcmp(all_intersection_high{2*k,i},'UU')
            NCloudU=length(all_intersection_high{2*k-1,i});
            index_1=i; 
            break;
        end
    end
       
    for i=1:length(all_intersection_low(k,:))
        if strcmp(all_intersection_low{2*k,i},'LL')
            NCloudL=length(all_intersection_low{2*k-1,i});
            index_2=i;
            break;
        end
    end   
 
    if(NCloudU==0||NCloudL==0)
        errordlg('Empty region');
        break;
    end

    Ncloud=min([NCloudU NCloudL]);
        
    MIRI=[];
     
    for ip=1:size(p_low,2) %for each parameter
        
        % parameter axis definition
        max_p=max(max(p_low(ip).*Results{k,1}(:,ip)),max(p_high(ip).*Results{k,1}(:,ip)));
        min_p=min(min(p_low(ip).*Results{k,1}(:,ip)),min(p_high(ip).*Results{k,1}(:,ip)));
        
        parameter_axis=(min_p:(max_p - min_p)/Ncloud:max_p);
                
        % kernel density estimation
        [ks_UUU xbin_UUU]=ksdensity(p_high(ip).*Results{k,1}(all_intersection_high{2*k-1,index_1}(:),ip),parameter_axis);
        [ks_LLL xbin_LLL]=ksdensity(p_low(ip).*Results{k,1}(all_intersection_low{2*k-1,index_2},ip),parameter_axis);        
       
        % plot probability density functions 
        figure(ip)
        plot(xbin_UUU, ks_UUU,'r',xbin_LLL,ks_LLL,'b');
        hold on;
               
        % MIRI computation      
        s_Xi_ks=sum(abs(ks_UUU./sum(ks_UUU)-ks_LLL./sum(ks_LLL)));        
        MIRI=[MIRI,s_Xi_ks]; 

    end
   
    MIRIT=[MIRIT;MIRI]; 
    
end
  
figure 
boxplot(MIRIT);

