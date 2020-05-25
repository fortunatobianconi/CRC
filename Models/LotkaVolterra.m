% Equations of the Lotka-Volterra model taken from:
% Toni, Tina, et al. 
% "Approximate Bayesian computation scheme for parameter inference and model selection in dynamical systems." 
% Journal of the Royal Society Interface 6.31 (2009): 187-202.


function dx=LotkaVolterra(t,x,p)

dx=zeros(2,1);


%%-----------------------------------%%
%              Equations              %
%%-----------------------------------%%

dx(1)=p(1)*x(1)-x(1)*x(2);
dx(2)=p(2)*x(1)*x(2)-x(2);


end