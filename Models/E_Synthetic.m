
% ODE of the model taken from:
% Raue, Andreas, et al. "Addressing parameter identifiability by model-based experimentation."
% IET systems biology 5.2 (2011): 120-130.


function dx=E_Synthetic(t,x,L,parameter)

k=parameter;


dx = zeros(5,1);

%%----------------%%
%    Parameters    %
%%----------------%%

k1=k(1);
k2=k(2);
k3=k(3);


%%------------------------------%%
%   Protein concentrations       %
%%------------------------------%%

E=x(1);
E_active=x(2);
E_active_active=x(3);
S=x(4);
P=x(5);


%%-----------------------------------%%
%              Equations              %
%%-----------------------------------%%

d_E = - k1*E*L;
d_E_active = + k1*E*L - k2*E_active;
d_E_active_active = k2*E_active;
d_S = -k3*E_active_active*S;
d_P= k3*E_active_active*S;


%%-----------------------------------%%
%           Differential Term         %
%%-----------------------------------%%

dx(1)=d_E;
dx(2)=d_E_active;
dx(3)=d_E_active_active;
dx(4)=d_S;
dx(5)=d_P;

end

