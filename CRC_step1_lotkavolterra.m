clear all
close all

tic;

% load the model to calibrate
addpath('Models')
model=@LotkaVolterra;
 
% time simulation
ts=0;
tf=15;   
T=1/100;
time=[ts:T:tf];

% initial conditions
ini_val=[1 0.5];

% model parameters to perturb
% in the first iteration, perturbed model parameters are all set to one,
% while in the subsequent iterations their nominal value is equal to the
% mode of the probability density function computed in the previous
% iteration and then the parameters are perturbed around the mode vector.
pN=[1 1];

%load('mode_parameters_iteration1.mat')
%pN=mode_parameters_iteration1;

% Latin hypercube perturbation of the parameter space
LBpi=0.1; % Lower bound
UBpi=10;   % Upper bound

Nsample=10000; % number of samples for each parameters
Nr=1; % number of independent realizations

% number of observables
Nobs=2;

% cell array where final results will be stored
Results=cell(Nr,Nobs+1);

% Performs perturbed simulations and store data
poolobj=parpool(4);
for k=1:Nr
    k
    Nominal_X1=zeros(8,Nsample);
    Nominal_X2=zeros(8,Nsample);
    
    % Logarithmic LHS sampling
    lb=log(LBpi).*[1 1];
    ub=log(UBpi).*[1 1];
    lhs=lhsdesign(Nsample,2);
    PerturbationLog=bsxfun(@plus, bsxfun(@times, lhs,ub-lb), lb);
    Perturbation(:,:)=exp(PerturbationLog(:,:));

    % Linear LHS sampling 
    % Perturbation=LBpi+(UBpi-LBpi).*lhsdesign(Nsample,2);
    
    Results{k,1}=Perturbation;
    
    parfor i=1:Nsample
        
        % model simulation for each generated parameter sample
        pi=pN.*Perturbation(i,:);
        [tsimi psimi]=ode15s(@(t,y) LotkaVolterra(t,y,pi), time, ini_val);
        msgstr=lastwarn;
        
        % storing of results
        if isequal(msgstr,'')
            i_X1=[psimi(111,1); psimi(241,1); psimi(391,1); psimi(561,1); psimi(751,1); psimi(961,1); psimi(1191,1); psimi(1441,1)];
            i_X2=[psimi(111,2); psimi(241,2); psimi(391,2); psimi(561,2); psimi(751,2); psimi(961,2); psimi(1191,2); psimi(1441,2)];
            
            Nominal_X1(:,i)=i_X1;
            Nominal_X2(:,i)=i_X2;
             
        else
            
            Nominal_X1(:,i)=[NaN; NaN; NaN; NaN; NaN; NaN; NaN; NaN];
            Nominal_X2(:,i)=[NaN; NaN; NaN; NaN; NaN; NaN; NaN; NaN];
                        
        end
        lastwarn('')
    end
    Results{k,2}=Nominal_X1;
    Results{k,3}=Nominal_X2;
     
end
delete(poolobj);
toc;
save('lotkavolterra_measures.mat');

% load dataset for computation of the distance functions
addpath('Data')
load('Data/data_LotkaVolterra.mat');

% computation of distance functions
distance_iteration1=cell(Nr,Nobs);
data=data_LotkaVolterra;
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

save('distance_iteration1.mat') 