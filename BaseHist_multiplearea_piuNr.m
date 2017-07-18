clear all
close all

% load measures and distance functions 
load('Measures.mat');
load('distance_iteration1.mat');

% observables names
proteinName={'Y1','Y2'};


% tolerances definition for all observables   
lowThr=[51 31];
highThr=[];

accepted_dist=upperLowerSet_Nr(distance_iteration1,Nr,Nsample,proteinName,lowThr,highThr)

% perform intersection among selected values in order to define the region
% of interest
all_region={};
 for k=1:Nr
  region=intersection(accepted_dist(2*k-1,:));
  all_region=[all_region;region];
 end

% store parameter density modes 
CloudCondLT=[];

% conditional parameter density estimation
xbin_p={};
ks_p={};

for k=1:Nr
    
    % selection of the desired region
    NCloud=0;
    for i=1:length(all_region(k,:))
        if strcmp(all_region{2*k,i},'LL')
            NCloud=length(all_region{2*k-1,i});
            index_2=i;
            break;
        end
    end   
        
    CloudCondL=[];   
       
    % definition of all the parameters to evaluate (e.g. kinetics, initial
    % conditions and scale factors)
    
    pN_all=[pN ini_val(index_non_zero_ini_val) scale_factors];
    
    % for each parameter       
    for ip=1:length(pN_all)
        
        Parameter_axis=min(pN_all(ip).*Results{k,1}(:,ip)):(max(pN_all(ip).*Results{k,1}(:,ip))-min(pN_all(ip).*Results{k,1}(:,ip)))/NCloud:max(pN_all(ip).*Results{k,1}(:,ip));
                
        % conditional parameter density estimation
        [ks_L, xbin_L]=ksdensity(pN_all(ip).*Results{k,1}(all_region{2*k-1,index_2},ip),Parameter_axis);        
        
        % store results
        xbin_p{k,ip}=xbin_L;
        ks_p{k,ip}=ks_L;
        
        % mode of the parameter density
        [ks_L_m, ks_L_I]=max(ks_L);
        
        CloudCondL=[CloudCondL,xbin_L(ks_L_I)];
        
     end
     CloudCondLT=[CloudCondLT;CloudCondL];
 end
  
