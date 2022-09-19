function out = CovFit2(c,t)
% ---------------------------------------------------------------input data
% (c,t) 
% length(t) x 3 data matrix. susceptible, infected, deceased


% DONT FORGET TO CHANGE RESULTS DIRECTORIES BEFORE EACH RUN
%      AND INITIAL CONDITIONS IN func1, func2, func4 etc.  
% ---------------------------------------------------------------

% gibb=fopen('C:\Users\pater\Documents\Spring 2021\COVID\Results 2\COVID\gibbs.dat','wt');

% rate guess
r0=[1;1;1;1];

[r]=lsqcurvefit(@func1,r0,t,c);
fprintf(1,'\tRate Constants:\n')

for k1 = 1:length(r)
    fprintf(1, '\t\tRate(%d) = %8.5f\n', k1, r(k1))
end

out='All Done';


function Conc=func2(r,t)
% initial conditions
C0=[1 0 0];
options = odeset('RelTol',1e-3);
[T,Concv]=ode45(@DifEq,t,C0,options);

function dConc=DifEq(t,C)

dCdt=zeros(4,1);

% system of equations
% susceptible
dCdt(1)=-C(1).*((r(1).*C(2)+r(3).*C(3)));
% infected A
dCdt(2)=r(1).*(C(1).*C(2))-r(2).*C(2);
% infected B
dCdt(3)=r(3).*(C(1).*C(3))-r(4).*C(3);
% dead
dCdt(4)=r(2).*C(2)+r(4).*C(3);


dConc=dCdt;

end
Conc=Concv;
end


end
