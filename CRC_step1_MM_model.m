clear all
close all

tic;

% load the model to calibrate
addpath('Models')
model=@MM_model; 

% time points of the simulation
time_points=[0 5 10 30 60 90];

% initial conditions
ini_val=ones(1,40); 

% definition of the model parameters and the corresponding non-zero
% elements.
% in the first iteration, perturbed model parameters are all set to one,
% while in the subsequent iterations their nominal value is equal to the
% mode of the probability density function computed in the previous
% iteration and then the parameters are perturbed around the mode vector.
pN=ones(1,53);
%load('mode_parameters_iterations1.mat');
%pN=mode_parameters_iterations1;

% array of the observables indexes in the experimental dataset
load('Data/exp_data_index.mat');
% array of the observables indexes in the ODE model 
load('Data/sim_data_index.mat');

% Latin hypercube perturbation of the parameter space
LBpi=0.01; % Lower bound
UBpi=100;   % Upper bound

Nsample=1000000; % number of samples for each parameters
Nr=10; % number of independent realizations

% number of observables
Nobs=16;

% cell array where final results will be stored
Results=cell(Nr,Nobs+1);

% Performs parallel perturbed simulations and store data
poolobj=parpool(4);
for k=1:Nr
    k
    
    % Logarithmic LHS sampling 
    lb=log(LBpi*ones(1,length(pN)));
    ub=log(UBpi*ones(1,length(pN)));
    lhs=lhsdesign(Nsample, length(pN));
    PerturbationLog=bsxfun(@plus,bsxfun(@times,lhs,ub-lb),lb);
    Perturbation(:,:)=exp(PerturbationLog(:,:));
    
    % Linear LHS sampling
    %Perturbation=LBpi+(UBpi-LBpi).*lhsdesign(Nsample,length(pN));
    
       
    Nominal_Shc=zeros(6,Nsample);
    Nominal_Pras=zeros(6,Nsample);
    Nominal_PI3K=zeros(6,Nsample);
    Nominal_p38=zeros(6,Nsample);
    Nominal_PDK1=zeros(6,Nsample);
    Nominal_Akt=zeros(6,Nsample);
    Nominal_mTOR=zeros(6,Nsample);
    Nominal_Raf=zeros(6,Nsample);
    Nominal_MEK=zeros(6,Nsample);
    Nominal_MAPK=zeros(6,Nsample);
    Nominal_p70=zeros(6,Nsample);
    Nominal_CJUN=zeros(6,Nsample);
    Nominal_Bclxl=zeros(6,Nsample);
    Nominal_BAX=zeros(6,Nsample);
    Nominal_NFKB=zeros(6,Nsample);
    Nominal_PARP=zeros(6,Nsample);
    
 
    Results{k,1}=Perturbation;

    parfor i=1:Nsample
        
       % model simulation for each generated parameter sample
       pi=pN.*Perturbation(i,:);
       
       [T y]=ode15s(@(t,y)ode_full(t,y,pi),time_points,ini_val);
        
       % definition of the output variables and results storage
        msgstr=lastwarn;
        if isequal(msgstr,'')
            sim_data_i=y(:,sim_data_index);
            
            Nominal_Shc(:,i)=sim_data_i(:,1);
            Nominal_Pras(:,i)=sim_data_i(:,2);
            Nominal_PI3K(:,i)=sim_data_i(:,3);
            Nominal_p38(:,i)=sim_data_i(:,4);
            Nominal_PDK1(:,i)=sim_data_i(:,5);
            Nominal_Akt(:,i)=sim_data_i(:,6);
            Nominal_mTOR(:,i)=sim_data_i(:,7);
            Nominal_Raf(:,i)=sim_data_i(:,8);
            Nominal_MEK(:,i)=sim_data_i(:,9);
            Nominal_MAPK(:,i)=sim_data_i(:,10);
            Nominal_p70(:,i)=sim_data_i(:,11);
            Nominal_CJUN(:,i)=sim_data_i(:,12);
            Nominal_Bclxl(:,i)=sim_data_i(:,13);
            Nominal_BAX(:,i)=sim_data_i(:,14);
            Nominal_NFKB(:,i)=sim_data_i(:,15);
            Nominal_PARP(:,i)=sim_data_i(:,16);          
            
            
        else
            
            Nominal_Shc(:,i)=nan(6,1);
            Nominal_Pras(:,i)=nan(6,1);
            Nominal_PI3K(:,i)=nan(6,1);
            Nominal_p38(:,i)=nan(6,1);
            Nominal_PDK1(:,i)=nan(6,1);
            Nominal_Akt(:,i)=nan(6,1);
            Nominal_mTOR(:,i)=nan(6,1);
            Nominal_Raf(:,i)=nan(6,1);
            Nominal_MEK(:,i)=nan(6,1);
            Nominal_MAPK(:,i)=nan(6,1);
            Nominal_p70(:,i)=nan(6,1);
            Nominal_CJUN(:,i)=nan(6,1);
            Nominal_Bclxl(:,i)=nan(6,1);
            Nominal_BAX(:,i)=nan(6,1);
            Nominal_NFKB(:,i)=nan(6,1);
            Nominal_PARP(:,i)=nan(6,1);                         
            
        end
        lastwarn('')
    end
    Results{k,2}=Nominal_Shc;
    Results{k,3}=Nominal_Pras;
    Results{k,4}=Nominal_PI3K;
    Results{k,5}=Nominal_p38;
    Results{k,6}=Nominal_PDK1;
    Results{k,7}=Nominal_Akt;
    Results{k,8}=Nominal_mTOR;
    Results{k,9}=Nominal_Raf;
    Results{k,10}=Nominal_MEK;
    Results{k,11}=Nominal_MAPK;
    Results{k,12}=Nominal_p70;
    Results{k,13}=Nominal_CJUN;
    Results{k,14}=Nominal_Bclxl;
    Results{k,15}=Nominal_BAX;
    Results{k,16}=Nominal_NFKB;
    Results{k,17}=Nominal_PARP;
     
end

delete(poolobj);

toc;

save ('measures_MM_model.mat','-v7.3');

% load dataset for computation of the distance functions
addpath('Data')
load('Data/data_MM_model.mat');

% computation of distance functions
distance_iteration1=cell(Nr,Nobs);
data=data_MM_model;
for k=1:Nr
    for i=1:Nobs
        for j=1:Nsample
             
             % NNDF
             distance_iteration1{k,i}(1,j)= (nansum(abs(data(:,i) - Results{k,i+1}(:,j))));
             
             % NDF 
             %distance_iteration1{k,i}(1,j)= 1/size(data,1) * sum(abs(data(:,i) - Results{k,i+1}(:,j))./data(:,i));
        end
    end
end
 
save('distance_iteration1_MM_model.mat') 