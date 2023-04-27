% Samoilov Artem, 21.11.2020
% Task 3 for System Analysis Practicum
% Unauthorizd use of creator's intellectual property is strictly forbidden

%% Subtask 1
clc
N = 50;
x = 1:N;
fn = @(n) ((-1).^n)./(n.^2);
Sn = zeros(1,N);
Psi = zeros(1,N);

Sn = fn(x(:));
Sn = cumsum(Sn);
Psi = 1./(x.^2);
S = (-1)*ones(1,N)*(pi^2/12);
plot(x, abs(S-Sn), x, Psi);
legend({'S-Sn','Psi'});
xlabel('x');
ylabel('y');

clear

%% Subtask 2
clc
a = -4;
b = 4;
n = 100;
x = linspace(a, b, n);
f = @(x) cos(x);
g = @(x) x./pi;
F = f(x);
G = g(x);

plot(x,F,x,G);
xlabel('x');
ylabel('y');
legend({'cos(x)','x/pi'});
hold on;

h = @(x) cos(x)-x./pi;
[X,Y] = ginput(1);
X1 = fzero(h, X);
plot(X1,f(X1),'*');
legend({'cos(x)','x/pi','zero'});
hold off;
X1
Residual = abs(h(X1))
clear

%% Subtask 3
clc
a = -2;
b = 2;
n = 10000;
x = linspace(a,b,n);
f = @(x) double(x ~= 0).*x.*sin(1./x);
y = f(x);
X = zeros(1,n);
for i = 1:n
    if(x(i) < -0.5)
        X(i) = -1/pi;
    elseif(x(i) > 0.5)
        X(i) = 1/pi;
    else
        X(i) = fzero(f,x(i));
    end
end
plot(x,X);
xlabel('x');
ylabel('roots');
clear


%% Subtask 4
clear
clc
t0 = 0;
t1 = 400;
alpha = 1;
e1 = [-0.4, 0.2];
v1 = [1, 1];
xy00 = [e1, v1];

eps = 0.0001;
options = odeset('Events', @boom, 'InitialStep', 0.00001);
t = t0;

x = linspace(-1,1,100);
y = sqrt(1-x.^2);
plot(x,y,'b', x,(-1)*y,'b', e1(1),e1(2), '*');
hold on
xlim([-1,1]);
ylim([-1,1]);
xlabel('x');
ylabel('y');
daspect([1 1 1]);
flag = 1;

while flag
    [tt, yy, te, ye, ie] = ode45(@f, [t0 t1], xy00, options);
    xy = [yy(:,1), yy(:,2)];
    comet(xy(:, 1), xy(:, 2));
    if((size(te))~=0)
        t0 = te;
    else
        t0 = t1;
        flag = 0;
    end
    if(size(ye,1)~=0)
        y = @(x) (ye(2)/ye(1))*x;
        xeps = ye(1) - (ye(1)/abs(ye(1)))*eps;
        yeps = y(xeps);
        start = [xeps, yeps];
        e1 = [ye(1), ye(2)];
        e2 = [-ye(2), ye(1)];
        v1 = [ye(3), ye(4)];
        v2 = e2.*dot(e2,v1)-(1/alpha).*e1.*dot(e1,v1);
        xy00 = [start, v2];
        hold on
    end
        
end
hold off

%% Subtask 5
clc
clear
T = [0, 1];
u0 = [100, 5];

[x, y] = ode45(@munch, T, u0);
sz = size(x,1);

for i = 1:sz
    plot(x(1:i),y(1:i,1), x(1:i),y(1:i,2));
    axis tight;
    legend({'prey','predators'});
    pause(0.001);
end
hold on

tt = linspace(0,1,10);
X = 100*exp(2*tt+exp(-10*tt)/2-1/2);
Y = 5*exp(-10*tt);
plot(tt,X,tt,Y);
legend({'prey','predators','an-prey','an-predators'});
hold off

%% Subtask 6
% узел
clc
clear
n = 50;
x = linspace(0, 2*pi, n);

a = [-3, 2; -2, 1];

g = @(t, xy) feval(@linsyst, t, xy, a);

T = [-10, 10];
r = 10;

for i = 1:n
    y0 = [r*cos(x(i)), r*sin(x(i))];
    [t, xy] = ode45(g, T, y0);
    plot(xy(:,1), xy(:,2));
    hold on;
end
x = linspace(-10, 10);
y = linspace(-15, 15);
[x, y] = meshgrid(x, y);
u = 2*y-3*x;
v = y - 2*x;
q = quiver(x, y, u, v);
daspect([1 1 1]);
hold off

%% Subtask 6
% центр
clc
clear
n = 50;
x = linspace(0, 2*pi, n);

a = [1, -1; 2, -1];

g = @(t, xy) feval(@linsyst, t, xy, a);

T = [-5, 5];
r = 4;

for i = 1:n
    y0 = [r*cos(x(i)), r*sin(x(i))];
    [t, xy] = ode45(g, T, y0);
    plot(xy(:,1), xy(:,2));
    hold on;
end
x = linspace(-10, 10);
y = linspace(-10, 10);
[x, y] = meshgrid(x, y);
u = x-y;
v = 2*x-y;
q = quiver(x, y, u, v);
daspect([1 1 1]);
hold off

%% Subtask 6
% седло
clc
clear
n = 50;
x = linspace(0, 2*pi, n);

a = [3, -4; 1, -2];

g = @(t, xy) feval(@linsyst, t, xy, a);

T = [-0.25, 0.25];
r = 1;

for i = 1:n
    y0 = [r*cos(x(i)), r*sin(x(i))];
    [t, xy] = ode45(g, T, y0);
    plot(xy(:,1), xy(:,2));
    hold on;
end
x = linspace(-5, 5);
y = linspace(-2, 2);
[x, y] = meshgrid(x, y);
u = 3*x-4*y;
v = x-2*y;
q = quiver(x, y, u, v);
daspect([1 1 1]);
hold off

%% Subtask 6
% устойчивый фокус
clc
clear
n = 50;
x = linspace(0, 2*pi, n);

a = [0, 1; -2, -1];

g = @(t, xy) feval(@linsyst, t, xy, a);

T = [-1, 1];
r = 8;

for i = 1:n
    y0 = [r*cos(x(i)), r*sin(x(i))];
    [t, xy] = ode45(g, T, y0);
    plot(xy(:,1), xy(:,2));
    hold on;
end
x = linspace(-10, 10);
y = linspace(-10, 10);
[x, y] = meshgrid(x, y);
u = y;
v = -y-2*x;
q = quiver(x, y, u, v);
daspect([1 1 1]);
hold off

%% Subtask 6
% неустойчивый фокус
clc
clear
n = 50;
x = linspace(0, 2*pi, n);

a = [0, 1; -2, 1];

g = @(t, xy) feval(@linsyst, t, xy, a);

T = [-0.2, 0.2];
r = 6;

for i = 1:n
    y0 = [r*cos(x(i)), r*sin(x(i))];
    [t, xy] = ode45(g, T, y0);
    plot(xy(:,1), xy(:,2));
    hold on;
end
x = linspace(-10, 10);
y = linspace(-10, 10);
[x, y] = meshgrid(x, y);
u = y;
v = -2*x+y;
q = quiver(x, y, u, v);
daspect([1 1 1]);
hold off

%% Subtask 7-1
clc
clear
v = @(x,y) x.^2+y.^2;
%dv/dt = -2(x-y)^4-2y^4
t0 = -10;
t1 = 10;
n = 20;
xx = linspace(0, 2*pi, n);
edge = 5;
[X, Y] = meshgrid(-edge:0.02:edge, -edge:0.02:edge);
Z = v(X, Y);
contour(X, Y, Z, 20);
hold on


r = 5;
for i=1:n
    y0 = [r*cos(xx(i)), r*sin(xx(i))];
    [t, F] = ode45(@system1, [t0, t1], y0);
    for j = 1:size(F,1)-1
        lin = plot(F(j:j+1, 1), F(j:j+1, 2),'LineWidth',2);
        lin.Color = [1 1/3+j/((3/2)*size(F,1)-1) 0];
        hold on
    end
end

hold on
x = linspace(-5, 5);
y = linspace(-5, 5);
[x, y] = meshgrid(x, y);
u = y-x+x.*y;
v = x-y-x.^2-y.^3;
q = quiver(x,y, u,v);
daspect([1 1 1]);

title('System 1: x''=y-x+xy; y''=x-y-x^2-y^3');
legend('Trajectory');
xlabel('x');
ylabel('y');
hold off


%% Subtask 7-2
clc
clear
v = @(x,y) x.^2+y.^2;
%x^2*y-2*y^4+x^2y^2
%dv/dt = 2x^3-2xy^3+2x^2y^2
%(x-y).^4;%+y.^2;
%v = @(x,y) x.*y.^3;
t0 = -1;
t1 = 1;
n = 50;

xx = linspace(0, 2*pi, n);
[X, Y] = meshgrid(-0.4:0.02:2, -0.4:0.02:0.4);
Z = v(X, Y);
contour(X, Y, Z, 20);
hold on

r = 0.4;
for i=1:n
    y0 = [r*cos(xx(i)), r*sin(xx(i))];
    [t, F] = ode45(@system2, [t0, t1], y0);
    for j = 1:size(F,1)-1
        lin = plot(F(j:j+1, 1), F(j:j+1, 2),'LineWidth',2); 
        lin.Color = [1 1-j/(size(F,1)-1) 0];
        hold on
    end
end
[x, y] = meshgrid(-0.4:0.02:2, -0.4:0.02:0.4);
u = x.^2-2*y.^3;
v = x.*y.^2;
q = quiver(x,y, u,v);
daspect([1 1 1]);

title('System 2: x''=x^2-2y^3; y''=xy^2');
legend('Trajectory');
xlabel('x');
ylabel('y');
hold off


%% Subtask 8
clc 
clear
x = linspace(0,pi/2,10);
solinit = bvpinit(x, @guess);
sol = bvp4c(@bvpfcn, @bcfcn, solinit);
sol.y(1,:);
analit = -cos(sol.x)-sin(sol.x)+1;
plot(sol.x, sol.y(1,:), sol.x, analit);
diffl2 = sqrt(sum((sol.y(1,:)-analit).^2))
diffC = abs(max(abs(sol.y(1,:)))-max(abs(analit)))

%% Subtask 9
clear
clc
f = @(x,y) x.^2.*y.^2+x.^2+y.^2;
gradx = @(x,y) 2*x.*y.^2+2*x;
grady = @(x,y) 2*x.*y.^2+2*y;
d = 0.0001;
x1 = 2;
y1 = 2;
k = 1;
kmax = 100000;
xtrace = [x1];
ytrace = [y1];
i = 2;
while k < kmax
    [xy, grad] = f2minbnd(gradx,grady,x1,y1);
    xtrace(i) = xy(1); 
    ytrace(i) = xy(2); 
    i = i + 1;
    x1 = xy(1);
    y1 = xy(2);
    if abs(grad) <= d
        break;
    end
    k = k + 1;
end
xtrace;
ytrace;
x = -2.5:0.01:2.5;
y = -2.5:0.01:2.5;
[X, Y] = meshgrid(x, y);
F = f(X,Y);

for i = 1:size(xtrace,2)
    Z = [f(xtrace(i),ytrace(i)), f(xtrace(i),ytrace(i))];
    [C, h] = contour(X,Y,F,Z);
    clabel(C, h);
    hold on;
end

plot(xtrace, ytrace, '-x');
text(xtrace(1), ytrace(1), 'XY1');
text(1, -2, char(['xf = ' (num2str(x1))], ['yf = ' (num2str(y1))], ['k = '  (num2str(k))]));
    
%% Subtask 9 - testing
clc
clear
fun = @(x) x.^2;
grad = @(x) 2*x;
x1 = 2;
g = 0.05;
xtrace = [x1];
i = 2;
d = 0.0001;
k = 1;
kmax = 10000;
while k < kmax 
    gr = grad(x1); 
    x1 = x1 - g*abs(gr);
    xtrace(i) = x1;
    i = i + 1; 
    if abs(gr) <= d 
        break;
    end
    k = k + 1; 
end
xf = fminbnd(@(x) x.^2, 0,2)
yf = fun(xf)
pl = plot(xf,yf,'*');
pl.Color = 'g';
hold on
plot(xtrace, fun(xtrace), '-x');
xtrace(end)
fun(xtrace(end))

%% Subtask 10-1
clc
clear

hFigure = figure();
subplot(2,1,1);
xlabel("Re(w)");
ylabel("Re(w)");
hold on
subplot(2,1,2);
xlabel("Im(w)");
ylabel("Im(F(w))");

step = 0.1;
inpLimVec = [-0.3 0.4];
outLimVec = [0 700];
%outLimVec = [-100*pi 100*pi];
S = plotFT(hFigure, func1, ftfunc1, step, inpLimVec, outLimVec)

%% Subtask 10-2
clc
clear

hFigure = figure();
subplot(2,1,1);
xlabel("Re(w)");
ylabel("Re(F(w))");
hold on
subplot(2,1,2);
xlabel("Im(w)");

ylabel("Im(F(w))");

step = 0.05;
inpLimVec = [-0.4 0.5];
outLimVec = [-15 15];
S = plotFT(hFigure, func2, ftfunc2, step, inpLimVec, outLimVec)

%% Subtask 10-3
clc
clear

hFigure = figure();
subplot(2,1,1);
xlabel("Re(w)");
ylabel("Re(F(w))");
hold on
subplot(2,1,2);
xlabel("Im(w)");
ylabel("Im(F(w))");
step = 0.1;
inpLimVec = [0.8 1.5];
outLimVec = [-100 100];
S = plotFT(hFigure, func3, [], step, inpLimVec, outLimVec)

%% Subtask 10-4
clc
clear

hFigure = figure();
subplot(2,1,1);
xlabel("Re(w)");
ylabel("Re(F(w))");
hold on
subplot(2,1,2);
xlabel("Im(w)");
ylabel("Im(F(w))");
step = 0.2;
inpLimVec = [-1 1];
outLimVec = [-5 5];
S = plotFT(hFigure, func4, [], step, inpLimVec, outLimVec)
%% Subtask 10-5
clc
clear

hFigure = figure();
subplot(2,1,1);
xlabel("Re(w)");
ylabel("Re(F(w))");
hold on
subplot(2,1,2);
xlabel("Im(w)");
ylabel("Im(F(w))");
step = 0.01;
inpLimVec = [-500 500];
outLimVec = [-10 10];
S = plotFT(hFigure, func5, [], step, inpLimVec, outLimVec)
%%
clc
clear

hFigure = figure();
subplot(2,1,1);
xlabel("Re(w)");
ylabel("Re(F(w))");
hold on
subplot(2,1,2);
xlabel("Im(w)");
ylabel("Im(F(w))");
step = 0.01;
inpLimVec = [-10 10];
outLimVec = [-15 15];
overlay(hFigure, func4, [], step, inpLimVec, outLimVec)
%%
x = -10:0.1:10
y = func2(x)
plot(x,y(x))
syms t w
f1 = (exp(-2*abs(t)))./(1 + (cos(t)).^2);
f = -(cos(4*t)^2 - 1)/(3*t^2);
fourier(f,t,w)