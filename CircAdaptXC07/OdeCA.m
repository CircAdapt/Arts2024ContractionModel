function FI=OdeCA(DfDtFu,dt,T,fp)
% function [FI]=OdeCA(DfDtFu,dt,TMax,fp)
% Solving system of ordinary differential equation with given steps,
% 2nd order, 1 calculation of derivative per time step
% Especially designed for CircAdapt
% 1st state variable fi(1) is time with [dfi/dt(1)=DfDtFu(1)]=1
% DfDtFu= string pointing to function df/dt
% dt= integration step
% T = duration of time integration interval
% OUTPUT: t=column vector of time points
% FI: matrix, rows -> time : columns -> state variables
% Theo Arts, Oct 1, 2020
dfdt   = str2func(DfDtFu); % setting derivative function
nt     = round(T/dt)+1;    % max number of time points
np     = numel(fp);        % number of state variables
fi     = zeros(nt,np);     % reserve memory space for output
h      = fp;               % h=output, set to initial condition
fi(1,:)= h ;               % copy to output matrix
df1    = dfdt(h')';        % derivative at 1st sample point
f      = h + df1*dt;       % 1st estimate next sample point
df     = dfdt(f')';        % derivative at 2nd sample point
g      = 0.5*(df+df1);     % 1st estimate 2nd sample point
h      = h + g*dt;         % value of state variable, 2nd order approx
fi(2,:)= h;                % copy to output matrix
it     = 3;                % Regular integration starts with 3rd point
while it<=nt               
    f = h + g*dt;          % 1st estimate of next sample point
    df= dfdt(f')';         % derivative at next sample point
    g = (g+2*df)/3;        % updated derivative at this sample pont
    h = h + g*dt;          % value of state variable, 2nd order approx
    fi(it,:)=h;            % copy to output matrix
    it=it+1;
end
FI=fi(1:nt,:);
end

