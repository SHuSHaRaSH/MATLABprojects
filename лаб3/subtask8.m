%% subtask8

xmesh = linspace(0,1);
answ = @(x) exp(1)/(1-exp(2))*(exp(-x)-exp(x))-2*x;
 
solinit = bvpinit(xmesh, @guess); 
sol = bvp4c(@bvpfcn, @bcfcn, solinit);
plot(sol.x, sol.y(1,:));
xlabel('x');
ylabel('y');
legend('solution');
 
normL2 = sqrt(trapz(xmesh, (answ(xmesh) - sol.y(1, :)).^2))
normC = max(abs(answ(xmesh) - sol.y(1, :)))
 
function dydx = bvpfcn(x, y)
    dydx = [y(2); 
            2*x+y(1)];
end
 
function res = bcfcn(ya,yb)
    res = [ya(1); 
           yb(1)+1];
end
 
function g = guess(x)
    g = [sin(x); cos(x)];
end