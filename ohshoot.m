clear all
clf

% Integrate the 2nd order ODE BVP by using the SHOOTFIRST method:

% We write SHOOTFIRST in this case such that each psi(1) and psi(2) have
% one unknown boundary. The reason why either psi(1) or psi(2) does not
% have *both* boundary conditions from the start is because one of the
% original boundary conditions of the 2nd-order original equation is a
% derivative condition (ie, a gradient bound). This is DIFFERENT from the
% setup used in "finitedifference.m."

global Da % Damkohler number
Da = 12.0;
yrange = [0 1]; % symmetric function
psiTrialRange = [0 1];
[y,psi] = shootfirst(@phif,@(idk1)[idk1,0],@(idk2)[1,idk2],.001,psiTrialRange,yrange,100);


function dpsidy = phif(y,psi)
global Da
dpsidy(1) = psi(2);
dpsidy(2) = psi(1)*Da;
dpsidy = dpsidy';
end

function [y,psi] = shootfirst(dpsidy,psi0,psi1,TOL,psiTrialRange,yrange,N)
didk = (psiTrialRange(2)-psiTrialRange(1))/N;
idk1 = psiTrialRange(1); y = linspace(yrange(1),yrange(2),N);
f = ode45(@phif,[0,1],idk1); % "psi0(idk)" is equivalent to "[idk 0]" here.
if abs(f(idk1)) > TOL
    idk1 = idk1-f(idk1)/(f(idk1)-f(idk1-didk))*didk; % didk is step size... how to implement adaptive step in this kind of algorithm? what is err?
    while abs(f(idk1)) < TOL
    idk1 = idk1 - f(idk1)/(f(idk1+didk)-f(idk1-didk))*2*didk;
    end
end
%{
idk2 = psiTrialRange(1);
f = ode45(@phif,[0 1],psi1(idk2)); % "psi1(idk)" is equivalent to "[idk 0]" here.
if abs(f(idk2)) > TOL
    idk2 = idk2-f(idk2)/(f(idk2)-f(idk2-didk))*didk; % didk is step size... how to implement adaptive step in this kind of algorithm? what is err?
    while abs(f(idk2)) < TOL
    idk2 = idk2 - f(idk2)/(f(idk2+didk)-f(idk2-didk))*2*didk;
    end
end
psi = -ode45(@phif,[0 1],psi0(idk1))*Da
%}
% Above code block is unnecessary as, in this use case, psi(1) = psi!!!
end