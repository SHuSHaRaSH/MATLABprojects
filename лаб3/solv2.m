
function  sol = solv2; 
    xmesh = linspace(0,pi/2,5);
    solinit = bvpinit(xmesh, @guess);
    sol = bvp4c(@bvpfcn, @bcfcn, solinit);
    disp(sol.y);
end




function dydx = bvpfcn(x,y) % equation to solve
dydx = zeros(2,1);
dydx = [y(2)
       -y(1)];
end
%--------------------------------
function res = bcfcn(ya,yb) % boundary conditions
res = [ya(1)
       yb(1)-2];
end
%--------------------------------
function g = guess(x) % initial guess for y and y'
g = [sin(x)
     cos(x)];
end