% catdisk_implicit.m by Nate Lynd
% This code will solve for the non-steady state reaction-diffusion in a
% catalytic disk problem that we've been using. 
clear all
% 1st: Let's start with a collocation in space... it's y in the notes but
% we'll call it x here:
Ny = 100; % This is the number of points
y = linspace(-1,1,Ny);
% The spacing between points
delta_y = y(3) - y(2);
% Damkoehler number
Da = 1;
 
% Time... we'll continually overwrite the last psi^n with psi^n+1 at all
% points j
delta_t = 0.4;
% For implicit, this (del) will affect accuracy, but not stability:
global del
del = delta_t/(delta_y)^2; 
disp(['delta_t/(delta_y)^2 = ' num2str(del) ])
% Initialize the initial condition at t=0:
psi_old = zeros(Ny-2,1);% Ny-2 because we don't need to calculate our 
                        % our boundary conditions at j=1, and Ny
psi_new = zeros(Ny-2,1);
% Boundary conditions affect the j=2, and j=Ny-1 points in our system
% The indexing is strange here because we lost j=1, and j=Ny in the
% original collocation due to boundary conditions. The solution is stored
% in an array that is Ny-2 long due to these two element dropping out of
% the problem.
psi_old(1) = 0 - del;   % first element
psi_old(Ny-2) = 0 - del;% last element
% start at 
t = 0.;
% Set up the tri-diagonals:
a = ones(Ny-2,1).*del;
b = ones(Ny-2,1).*(-2*del-1-Da*delta_t);
d = ones(Ny-2,1).*del;
 
hold on
for tcalc = 0.0:0.1:1.0 % calculate solution at several t values
    while t < tcalc % integrate forward until t ~ tcalc, then add to plot
    % No need to update a,b, or d inside loop
        psi_new = tridiag(a,b,d,psi_old);
        t = t + delta_t;
        % get psi_old ready for the next iteration
        psi_old = -1.*psi_new;
        % reset the first and last elements
        psi_old(1) = -psi_new(1) - del;
        psi_old(Ny-2) = -psi_new(Ny-2) - del;
    end
    % appropriate place to add to a plot
    plot(y(2:Ny-1),psi_new)
end
hold off
legend('0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0')
xlabel('y - non-dimensionalized coordinate across the catalytic disk/plate')
ylabel('psi - non-dimensionalized concentration with Da = 1')

function psi_new = tridiag(a,b,d,psi_old)
global del
psi_new = zeros(length(b),length(b));
for i=1:1:length(psi_old)
    psi_old(i+1).*del-b.*psi_old+psi(i-1).*del
    psi_new
    psi_new(i) = psi_old(i+1).*del-b.*psi_old+psi(i-1).*del;
    psi_new(i) = -psi_new(i);
end
end