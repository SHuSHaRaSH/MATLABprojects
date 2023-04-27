%% task3
% %% 1
% 
% 
% a = -2;
% b = 0;
% x = a:0.01:b;   
% f1 = @(x) sin(x);
% f2 = @(x) x.^3 + x + 1;
% f3 = @(x) x.^3 + x + 1-sin(x);
% plot(x, f1(x), x, f2(x));
% hold on
% while 1
%     
%     [x1,y1] = ginput(1);
%     disp("нач реш =");
%     disp(x1);
%     plot(x1,y1,'O');
%     y1 = fzero(f3, x1); 
%     disp("найд реш =");
%     disp(y1);
%     plot(y1,f1(y1),'*');
%     disp("невязка=");
%     disp(abs(f3(y1)));
% end

% %% 2
%  clear
%     clc
%      f1 = @f;
%     a = 1;
%     x = -a:0.01:a;                      
%     plot(x, f1(x));                    
%     hold on
%              
%     
%     y = zeros(1, size(x,2));    
%     for i = 1:size(x,2)
%         y(i) = fzero(f1, x(i));               
%     end;
%     plot(x,y);
%     hleg = legend('x*cos(ln(|x|))', 'начальное приближение', 'корень по оси абсцисс', 'корень по оси ординат');
%     hold off

% %% 3
% clc,clear
%     tspan = [0 1];
%     A = [3 1; 4 1]
%     H = [A zeros(2,2); zeros(2,2) A];
%     y01 = [1 ;0;0;1];
%     [t, y] = ode45(@(t,y) H*y, tspan, y01);
%     disp(y);
%     expA = [transpose(y(end,1:2)), transpose(y(end,3:4))];
%     sum = [1 0; 0 1];
%     n = 10; 
%     k =1;
%     for i = 1:n
%         k = k*i;
%         sum =sum + (A^i)/k;
%     end
%     disp(expm(A));
%     disp(sum);
%     disp(expA);
%     
%     disp(abs(expA - expm(A)));
%     disp(abs(sum - expm(A)));
%     disp(abs(sum - expA));
    

   %%
   n = 15;
   x = 1:n;
   y = 1:n;
   
   mov(1:n) = struct('cdata', [], 'colormap', []);
   for i = 1:n
       hold on
       viscircles([x(i), y(i)], 1);
       hold off
       xlim([0 20]);
       ylim([0 20]);
       daspect([1 1 1]);
       mov(i) = getframe();
       clf
   end
   movie(mov)
   
   %% 4
   clc, clear;
   eps = 0.0001;
b = 10;
a = 5;
border = polyshape([-a a a -a], [-b -b b b]);
plot(border);
axis manual;
hold on;
time0 = 0;
time1 = 10;
init_cond = [1 -10 1 10];
s = (abs(init_cond(2))+abs(init_cond(4)))/100;
plot(init_cond(1), init_cond(3), 'ro');
options = odeset('Events', @eventfcn);
[t, x, te, xe, ie] = ode45(@odefcn, [time0:s:time1], init_cond, options);

comet(x(:, 1), x(:, 3));
disp(x(:, 3));
disp(x(:, 1));
disp(te);
disp(xe);
disp(ie);
while t(end) < time1
if ie == 1 && xe(1) < 0
init_cond = [xe(1)+eps -xe(2) xe(3) xe(4)];
elseif ie == 1 && xe(1) > 0
init_cond = [xe(1)-eps -xe(2) xe(3) xe(4)];
elseif ie == 2 && xe(3) < 0
init_cond = [xe(1) xe(2) xe(3)+eps -xe(4)];
elseif ie == 2 && xe(3) > 0
init_cond = [xe(1) xe(2) xe(3)-eps -xe(4)];
end
[t, x, te, xe, ie] = ode45(@odefcn, [t(end):s  :time1], init_cond, options);
comet(x(:, 1), x(:, 3));
end

  %% 5
  
 %% 6
%% 7
%% 8
%clear;

answ = @(x) exp(1)/(1-exp(2))*(exp(-x)-exp(x))-2*x;
 
solinit = bvpinit(xmesh, @guess); 
sol = bvp4c(@bvpfcn, @bcfcn, solinit);
plot(sol.x, sol.y(1,:));
xlabel('x');
ylabel('y');
legend('sov');
 
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

   