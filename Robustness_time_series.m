
clear all
close all

tic;

% load the model to calibrate
addpath('Models')
model=@E_Synthetic; 

% time points of the simulation
time_points=[0 3.33 6.66 10 13.33 16.66 20 23.33 26.66 30 40];

% definition of the model parameters and the corresponding non-zero elements 
ini_val=[1 0 0 1 0];
index_non_zero_ini_val=[1 4];

pN=[1 1 1];

scale_factors=[1 1];

% definition of the fixed input model parameters
L=1;

% Latin hypercube perturbation of the parameter space
LBpi=0.01; % Lower bound
UBpi=100;   % Upper bound

Nsample=1000; % number of samples for each parameter
Nr=10; % number of independent realizations

% number of observables
Nobs=2;

% cell array where final results will be stored
Results=cell(Nr,Nobs+1);

% Perform parallel perturbed simulations and store data

matlabpool open local 4
for k=1:Nr
    k
    Nominal_Y1=zeros(length(time_points),Nsample);
    Nominal_Y2=zeros(length(time_points),Nsample);
    
    % Logarithmic LHS sampling
    lb_pN=log(LBpi).*ones(1,length(pN));
    ub_pN=log(UBpi).*ones(1,length(pN));
    lhs=lhsdesign(Nsample,length(pN));
    PerturbationLog=bsxfun(@plus, bsxfun(@times, lhs,ub_pN-lb_pN), lb_pN);
    Perturbation_pN(:,:)=exp(PerturbationLog(:,:));
    
    lb_ini_val=log(LBpi).*ones(1,sum(ini_val~=0));
    ub_ini_val=log(UBpi).*ones(1,sum(ini_val~=0));
    lhs=lhsdesign(Nsample, sum(ini_val~=0));
    PerturbationLog=bsxfun(@plus, bsxfun(@times, lhs,ub_ini_val-lb_ini_val), lb_ini_val);
    Perturbation_ini_val(:,:)=exp(PerturbationLog(:,:));
    
    lb_scale_factors=log(LBpi).*ones(1,length(scale_factors));
    ub_scale_factors=log(UBpi).*ones(1,length(scale_factors));
    lhs=lhsdesign(Nsample,length(scale_factors));
    PerturbationLog=bsxfun(@plus, bsxfun(@times, lhs,ub_scale_factors-lb_scale_factors), lb_scale_factors);
    Perturbation_scale_factors(:,:)=exp(PerturbationLog(:,:));
    
    % Linear LHS sampling
    %Perturbation=LBpi+(UBpi-LBpi).*lhsdesign(Nsample,length(pN));
    
   
    Results{k,1}=[Perturbation_pN Perturbation_ini_val Perturbation_scale_factors];
    
    parfor i=1:Nsample
        
        % model simulation for each parameter sample
        pi=pN.*Perturbation_pN(i,:);
        
        ini_val_temp=ini_val(index_non_zero_ini_val).*Perturbation_ini_val(i,:);
        ini_val_i=ini_val;
        for j=1:length(index_non_zero_ini_val)
            ini_val_i(index_non_zero_ini_val(j))=ini_val_temp(j);
        end    
        
        [tsimi psimi]=ode15s(@(t,y) model(t,y,L,pi), time_points, ini_val_i);
        msgstr=lastwarn;
        
        % definition of the output variables
        scale_factors_i=scale_factors.*Perturbation_scale_factors(i,:);
        y=zeros(length(time_points),Nobs);
        y(:,1)=scale_factors(1)*psimi(:,1);
        y(:,2)=scale_factors(2)*psimi(:,4);
        
        % storing of results
        if isequal(msgstr,'')
            i_Y1=y(:,1);
            i_Y2=y(:,2);          
            Nominal_Y1(:,i)=i_Y1;
            Nominal_Y2(:,i)=i_Y2;  
        else            
            Nominal_Y1(:,i)=[NaN; NaN; NaN; NaN; NaN; NaN; NaN; NaN; NaN; NaN; NaN];
            Nominal_Y2(:,i)=[NaN; NaN; NaN; NaN; NaN; NaN; NaN; NaN; NaN; NaN; NaN];                        
        end
        lastwarn('')
    end
    Results{k,2}=Nominal_Y1;
    Results{k,3}=Nominal_Y2;
     
 end
matlabpool close
toc;

save('Measures.mat');

% load dataset for computation of the distance functions
load('Data/data_E_Synthetic.mat');

% distance functions
distance_iteration1=cell(Nr,Nobs);
data=data_E_Synthetic;
for k=1:Nr
    for i=1:Nobs
        for j=1:Nsample
             
             % NNDF
             distance_iteration1{k,i}(1,j)= (nansum(abs(data(:,i) - Results{k,i+1}(:,j))));
             
             % NDF 
             %error_iteration1{k,i}(1,j)= 1/size(data,1) * sum(abs(data(:,i) - Results{k,i+1}(:,j))./data(:,i));
        end
    end
end
 
save('distance_iteration1.mat') 
