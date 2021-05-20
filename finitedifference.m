clear all
clf

% Integrate the 2nd order ODE BVP by using the FINITE DIFFERENCE method:
[x1, y1] = uncool(.2); [x2, y2] = uncool(2.0); [x3, y3] = uncool(12.0);
[x4, y4] = uncool(120.0);
hold on
plot(x1,y1,'-r',x2,y2,'-b',x3,y3,'-g',x4,y4,'-y');
hold off

function [psi,y] = uncool(Da)
psi0 = 1; psiNplus1 = 1; N = 200;
yrange=[-1 1];
y = linspace(yrange(1),yrange(2),N); dy = y(2)-y(1);
d(1) = -1*psi0; d(N) = -1*psiNplus1; % Matlab automatically fills with zeros.
d = d';
a = 1; b = -(dy^2*Da+2); c = 1; % Doesn't change with i. But, if they did, we would make fa -> fa(i).
psi = d\tridiag(a,b,c,N,N);
d
end

function mat = tridiag(a,b,c,ni,nj)
mat = zeros(ni,nj);
mat(1,1:2) = [b,c]'; % first row, where d(1) = -a(1)
mat(ni,nj-1:nj) = [a,b]'; % last row, where d(N) = -c(N)
for i=2:1:ni-1
jstart = i-1;
% Raster through return mat
mat(i,jstart:jstart+2) = [a b c];
end
end

function [y,psi] = finite_diff(d2sysdy2,N,psi0andNp1)
% Lesson opportunity: Writing contained methods for 2nd order DE's is
% exceedingly complicated compared to 1st order DE's. Don't do this for FD
% method.
y = 0; psi = 0;
end

