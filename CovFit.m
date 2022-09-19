function out = CovFit(c,t)
% ---------------------------------------------------------------input data
% (c,t) 
% length(t) x 3 data matrix. susceptible, infected, deceased


% DONT FORGET TO CHANGE RESULTS DIRECTORIES BEFORE EACH RUN
%      AND INITIAL CONDITIONS IN func1, func2, func4 etc.  
% ---------------------------------------------------------------

% open results file for delta G
gibb=fopen('C:\Users\pater\Documents\Multistrain\gibbs.dat','wt');

% initial rate guess
r0=[2;1];

% curve fitting 
options = optimoptions(@lsqcurvefit,'OptimalityTolerance',1.0000e-6,'FunctionTolerance',1.0000e-6,'StepTolerance',1.0000e-6);
[r]=lsqcurvefit(@func1,r0,t,c,[],[],options);
fprintf(1,'\tRate Constants:\n')

for k1 = 1:length(r)
    fprintf(1, '\t\tRate(%d) = %8.5f\n', k1, r(k1))
end

t52=linspace(0,52,52);
Gibbs=CovG(rrow,c,tv);

for kk=1:52
fprintf(gibb, ' %3.5f \t %3.5f \t %3.5f \r\n', [t(kk); Gibbs(kk,1); Gibbs(kk,2)]);
end

function out=CovG(r,c,t)   
out=zeros(52,2);
for k=1:1:52

%forces
out(k,1)=-log(r(1).*c(k,1).^(c(k,2)/c(k,1)))+log(r(2).*c(k,2));
out(k,2)=-log(r(2).*c(k,2)^(c(k,3)/c(k,2)));
end

end


out='All Done';

function Conc=func1(r,t)

% susceptible, infected, deceased initial conditions
C0=[1,0,0];
options = odeset('RelTol',1e-6);
[T,Concv]=ode45(@DifEq,t,C0,options);

function dConc=DifEq(t,C)

dCdt=zeros(3,1);

% system of equations
% susceptible
dCdt(1)=-r(1).*(C(1).*C(2));
% infected
dCdt(2)=r(1).*(C(1).*C(2))-r(2).*C(2);
% deaceased
dCdt(3)=r(2).*C(2);

dConc=dCdt;

end
Conc=Concv;
end


end
