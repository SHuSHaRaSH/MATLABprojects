%%
clc
clear

%mex sbscale.c
%%
A = [1,1;1,1];
B = [1,2;1,3];
bias = 3;
scale = 2;
Y1 = sbscale(bias,scale,A,B);
Y2 = sbscalem(bias,scale,A,B);
Y1 == Y2

%% Subtask 1
clc
clear
% 1 -1.62 0.62
%A = complex([1],[1]);
%B = complex([-1],[-1]);
%C = complex([0],[0]);

%A = complex([-13],[-10]);
%B = complex([-19],[-1]);
%C = complex([-7],[-2]);

A = complex([1 2 19 1], [10 11 13 4]);
B = complex([4 -5 7 2], [13 14 0 5]);
C = complex([7 8 2 3], [16 17 -1 6]);

%A = complex([1, -2; 4, 8], [5, 11; 13, 5]);
%B = complex([-3, 10; 6, 5], [-1, 9; -7, 9]);
%C = complex([11, 12; 4, 1], [-5, -7; 9, 19]);

[X1, X2, X3] = cubesolvee(A,B,C);
X1
X2
X3
f = @(x) A.*(x.^3)+B.*x+C;
%f2 = @(x) A.*(imag(x).^3)+B.*imag(x)+C;
f(X1)
f(X2)
f(X3)

%% Subtask 2
f = @(x) x.*cos(x);
a = 0;
b = 2*pi;
x = linspace(a,b,30)';

[A1, B1, C1, D1] = createsplinee_c(x, f(x));
[A2, B2, C2, D2] = createsplinee_m(x, f(x));
cubespline = @(Ai, Bi, Ci, Di, x, xi) Ai + Bi*(x - xi) + 0.5*Ci*(x - xi).^2 + (1/6)*Di*(x - xi).^3;
hold on
plot(x, f(x));
for i = 2:size(x,1)
    linsp2 = linspace(x(i - 1), x(i), 10);
    plot(linsp2, cubespline(A1(i - 1), B1(i - 1), C1(i - 1), D1(i - 1), linsp2, x(i)), 'yellow');
    plot(linsp2, cubespline(A2(i - 1), B2(i - 1), C2(i - 1), D2(i - 1), linsp2, x(i)), 'red');
end
legend({'xcosx','.cpp','.m'});
xlabel('x');
ylabel('y');
hold off;

%% Subtask 3
f = @(x) x.*cos(x);

a = 0;
b = 2*pi;
x = linspace(a, b, 5);

linsp = [10, 20, 50, 70, 100, 150, 200, 300, 500, 700, 1000];
error = zeros(4, size(linsp, 2));

cubespline = @(Ai, Bi, Ci, Di, x, xi) Ai + Bi*(x - xi) + 0.5*Ci*(x - xi).^2 + (1/6)*Di*(x - xi).^3;

for i = 1 : size(linsp, 2)
    linsp2 = linspace(a, b, linsp(i));
    
    interp1f = interp1(x, f(x), linsp2, 'spline');
    error(1,i) = max(abs(f(linsp2)-interp1f));
    
    splinef = spline(x, f(x), linsp2);
    error(2,i) = max(abs(f(linsp2)-splinef));
    
    [A1, B1, C1, D1] = createsplinee_c(x',(f(x))');
    [A2, B2, C2, D2] = createsplinee_m(x',(f(x))');
    
    for j = 2:size(x,2)
        linspcurr = linsp2(linsp2 > x(j - 1) & linsp2 < x(j));
        errorcurr1 = max(abs(f(linspcurr) - cubespline(A1(j - 1), B1(j - 1), C1(j - 1), D1(j - 1), linspcurr, x(j))));
        errorcurr2 = max(abs(f(linspcurr) - cubespline(%%
clc
clear

%mex sbscale.c
%%
A = [1,1;1,1];
B = [1,2;1,3];
bias = 3;
scale = 2;
Y1 = sbscale(bias,scale,A,B);
Y2 = sbscalem(bias,scale,A,B);
Y1 == Y2

%% Subtask 1
clc
clear
% 1 -1.62 0.62
%A = complex([1],[1]);
%B = complex([-1],[-1]);
%C = complex([0],[0]);

%A = complex([-13],[-10]);
%B = complex([-19],[-1]);
%C = complex([-7],[-2]);

A = complex([1 2 19 1], [10 11 13 4]);
B = complex([4 -5 7 2], [13 14 0 5]);
C = complex([7 8 2 3], [16 17 -1 6]);

%A = complex([1, -2; 4, 8], [5, 11; 13, 5]);
%B = complex([-3, 10; 6, 5], [-1, 9; -7, 9]);
%C = complex([11, 12; 4, 1], [-5, -7; 9, 19]);

[X1, X2, X3] = cubesolvee(A,B,C);
X1
X2
X3
f = @(x) A.*(x.^3)+B.*x+C;
%f2 = @(x) A.*(imag(x).^3)+B.*imag(x)+C;
f(X1)
f(X2)
f(X3)

%% Subtask 2
f = @(x) x.*cos(x);
a = 0;
b = 2*pi;
x = linspace(a,b,30)';

[A1, B1, C1, D1] = createsplinee_c(x, f(x));
[A2, B2, C2, D2] = createsplinee_m(x, f(x));
cubespline = @(Ai, Bi, Ci, Di, x, xi) Ai + Bi*(x - xi) + 0.5*Ci*(x - xi).^2 + (1/6)*Di*(x - xi).^3;
hold on
plot(x, f(x));
for i = 2:size(x,1)
    linsp2 = linspace(x(i - 1), x(i), 10);
    plot(linsp2, cubespline(A1(i - 1), B1(i - 1), C1(i - 1), D1(i - 1), linsp2, x(i)), 'yellow');
    plot(linsp2, cubespline(A2(i - 1), B2(i - 1), C2(i - 1), D2(i - 1), linsp2, x(i)), 'red');
end
legend({'xcosx','.cpp','.m'});
xlabel('x');
ylabel('y');
hold off;

%% Subtask 3
f = @(x) x.*cos(x);

a = 0;
b = 2*pi;
x = linspace(a, b, 5);

linsp = [10, 20, 50, 70, 100, 150, 200, 300, 500, 700, 1000];
error = zeros(4, size(linsp, 2));

cubespline = @(Ai, Bi, Ci, Di, x, xi) Ai + Bi*(x - xi) + 0.5*Ci*(x - xi).^2 + (1/6)*Di*(x - xi).^3;

for i = 1 : size(linsp, 2)
    linsp2 = linspace(a, b, linsp(i));
    
    interp1f = interp1(x, f(x), linsp2, 'spline');
    error(1,i) = max(abs(f(linsp2)-interp1f));
    
    splinef = spline(x, f(x), linsp2);
    error(2,i) = max(abs(f(linsp2)-splinef));
    
    [A1, B1, C1, D1] = createsplinee_c(x',(f(x))');
    [A2, B2, C2, D2] = createsplinee_m(x',(f(x))');
    
    for j = 2:size(x,2)
        linspcurr = linsp2(linsp2 > x(j - 1) & linsp2 < x(j));
        errorcurr1 = max(abs(f(linspcurr) - cubespline(A1(j - 1), B1(j - 1), C1(j - 1), D1(j - 1), linspcurr, x(j))));
        errorcurr2 = max(abs(f(linspcurr) - cubespline(A2(j - 1), B2(j - 1), C2(j - 1), D2(j - 1), linspcurr, x(j))));
        
        if (j>2) 
            if errorcurr1 > error(3,i)
                error(3,i) = errorcurr1;
            end
            if errorcurr2 > error(4,i)
                error(4,i) = errorcurr2;
            end
        else
            error(3,i) = errorcurr1;
            error(4,i) = errorcurr2;
        end
    end
end

hold on;
plot(linsp, error(1, :));
plot(linsp, error(2, :));
plot(linsp, error(3, :));
plot(linsp, error(4, :));
legend({'interp1', 'spline', '.cpp', '.m'});
xlabel('number of points');
ylabel('error');
hold off;

%% Subtask 4

f = @(x) exp(x).*sin(x);

a = -pi;
b = pi;
x = linspace(a, b, 5);
linsp = 5:5:5000;
times = zeros(4, size(linsp,2));

repeats = 10;

cubespline = @(Ai, Bi, Ci, Di, x, xi) Ai + Bi*(x-xi) + 0.5*Ci*(x-xi).^2 + (1/6)*Di*(x-xi).^3;

for i = 1 : size(linsp,2)
    linsp2 = linspace(a,b,linsp(i));
    
    tic
    for k = 1:repeats
        interp1f = interp1(x,f(x),linsp2,'spline');
    end
    times(1,i) = toc;
    times(1,i) = times(1,i)/repeats;
    
    tic
    for k = 1:repeats 
        splinef = spline(x,f(x),linsp2);
    end
    times(2,i) = toc;
    times(2,i) = times(2,i)/repeats;
    
    tic
    for k = 1:repeats
        [A1, B1, C1, D1] = createsplinee_c(x',(f(x))');
    end
    times(3,i) = toc;
    times(3,i) = times(3,i)/repeats;
    
    tic
    for k = 1:repeats
        [A2, B2, C2, D2] = createsplinee_m(x',(f(x))');
    end
    times(4,i) = toc;
    times(4,i) = times(4,i)/repeats;
end

hold on;
plot(linsp, times(1,:));
plot(linsp, times(2,:));
plot(linsp, times(3,:));
plot(linsp, times(4,:));
legend({'interp1', 'spline', '.cpp', '.m'}, 'Location', 'northwest');
xlabel('number of points');
ylabel('execution time');
hold off;

%% Subtask 5
poldeg = 20;
polval = zeros(4, poldeg + 1);
polval(1,:) = polyfit(linsp, times(1,:), poldeg);
polval(2,:) = polyfit(linsp, times(2,:), poldeg);
polval(3,:) = polyfit(linsp, times(3,:), poldeg);
polval(4,:) = polyfit(linsp, times(4,:), poldeg);

for i = 1:4
    hold on;
    plot(linsp, times(i, :));
    plot(linsp, polyval(polval(i,:), linsp));
end
legend({'interp1', 'spline', '.cpp', '.m'}, 'Location', 'northwest');
xlabel('number of points');
ylabel('execution time');
hold off

%% Task 6-1
clc
clear
N = 100;
M = 100;
fHandle1 = fGiven;

mu = 1;
u1zero = 1;
u2zero = 1;
valnum = uNumerical(u1zero,u2zero,mu,N,M);
x = linspace(0,1-1/N,N);
y = linspace(0,1-1/M,M);
[X,Y] = meshgrid(y,x);
valanalit = uAnalytical(X,Y,u1zero,u2zero,mu);
figure;
hold on
surf(X,Y, valnum,'FaceAlpha',0.5,'FaceColor','g');
surf(X,Y, valanalit,'FaceAlpha',0.5,'FaceColor','b');
legend("Numerical","Analytical");
daspect([1 1 1]);
hold off
view(3);

figure;
hold on
surf(X,Y, abs(valanalit-valnum),'FaceAlpha',0.5,'FaceColor','y');
abs(valanalit-valnum)
%daspect([1 1 1]);
%legend("Numerical","Analytical","Error");
hold off
view(3);

%% Subtask 6-2
clc
clear
N = 100;
M = 100;
mu = 1;

% x(1-x)+y(1-y)
fHandle = @(x,y) -4 - mu.*(x.*(1-x)+y.*(1-y));
xiHandle = @(x) x.*(1-x);
etaHandle = @(y) y.*(1-y);

% x^2*sin(pi*x)+sin(4*pi*y)
%fHandle = @(x,y) (2-pi^2*x.^2).*sin(pi*x)+4*pi*x.*cos(pi*x)-16*pi.^2*sin(4*pi*y) - mu.*(x.^2.*sin(pi*x)+sin(4*pi*y));
%xiHandle = @(x) x.^2.*sin(pi*x);
%etaHandle = @(y) sin(4*pi*y);

% sinh(sin(pi*x))+sinh(sin(pi*y))
%fHandle = @(x,y) pi^2.*((cos(pi*x)).^2.*sinh(sin(pi*x))-sin(pi*x).*cosh(sin(pi*x)))+pi^2.*((cos(pi*y)).^2.*sinh(sin(pi*y))-sin(pi*y).*cosh(sin(pi*y)))- mu.*(sinh(sin(pi*x))+sinh(sin(pi*y)));
%xiHandle = @(x) sinh(sin(pi*x));
%etaHandle = @(y) sinh(sin(pi*y));

valnum = solveDirichlet2(fHandle,xiHandle,etaHandle,mu,N,M);
x = linspace(0,1-1/N,N);
y = linspace(0,1-1/M,M);
[X,Y] = meshgrid(y,x);
u =@(x,y) x.*(1-x)+y.*(1-y);
%u = @(x,y) x.^2.*sin(pi*x)+sin(4*pi*y);
%u = @(x,y) sinh(sin(pi*x))+sinh(sin(pi*y));

valanalit = u(X,Y);
hold on
surf(X,Y,real(valnum),'FaceAlpha',0.5,'FaceColor','g');
surf(X,Y,valanalit,'FaceAlpha',0.5,'FaceColor','b');
surf(X,Y,abs(valanalit-valnum),'FaceColor','y');

disp(max(max(abs(valnum-valanalit))))
legend("Numerical","Analytical","Error");
daspect([1 1 1]);
hold off
view(3);

%%
a = [1, NaN, NaN]
isNan(a)
sum(a)
