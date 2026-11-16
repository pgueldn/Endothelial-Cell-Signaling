function [tspan,y]= ODE_master(rpar,y0,speciesNames,tau,ymax, numEdges,ODEfilename,odeOptions)

tspan=linspace(0,35,100); %how long the ODE will run for to reach steady state 
 
odefun = @(t,y) ode_default(t,y,speciesNames, rpar,tau,ymax,numEdges,ODEfilename);
if nargin < 8 || isempty(odeOptions)
    [t,y] = ode45(odefun,tspan,y0);
else
    [t,y] = ode45(odefun,tspan,y0,odeOptions);
end
y=real(y); 
end

function dydt=ode_default(t,y,speciesNames,rpar,tau,ymax, numEdges,ODEfilename)
 
mapVar = cell(length(speciesNames),1);
for i = 1:length(speciesNames)
    mapVar{i} = sprintf('%s = %d;',speciesNames{i},i);
end
mapVar=horzcat(mapVar{1:end});
eval(mapVar); %Assign Stress=1 etc

load(ODEfilename);

ODElist=horzcat(ODElist{1:end});
eval(ODElist);

end

%% utility functions
function fact = act(x,rpar)
    w = rpar(1);
    n = rpar(2);
    EC50 = rpar(3); 
    beta = (EC50.^n - 1)./(2*EC50.^n - 1);
    K = (beta - 1).^(1./n);
 fact = min(w,w.*(beta.*x.^n)./(K.^n + x.^n));
end

function finhib = inhib(x,rpar)
    finhib = rpar(1) - act(x,rpar);
end

function z = OR(x,y)
    z = x + y - x*y;
end
    
function z = AND(rpar,varargin)
w=rpar(1);
   if mean(cell2mat(varargin))==0 % If all inputs are 0
            z=0;
         else
            v = cell2mat(varargin);
            z = prod(v)/w^(nargin-2); 
    end
end
