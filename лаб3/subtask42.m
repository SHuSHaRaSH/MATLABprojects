 eps = 0.0001;
b = 2;
a = 3;
border = polyshape([-a a a -a], [-b -b b b]);
plot(border);
axis manual;
hold on;
time0 = 0;
time1 = 30;
init_cond = [-1 5 0 10];
plot(init_cond(1), init_cond(3), 'ro');
options = odeset('Events', @eventfcn);
alpha = 0;
f = @(t, x) [x(3); x(4); -alpha*x(1); -alpha*x(2)];
[t, x, te, xe, ie] = ode45(f, [time0 time1], init_cond, options);
comet(x(:, 1), x(:, 3));
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
[t, x, te, xe, ie] = ode45(@odefcn, [t(end) time1], init_cond, options);
comet(x(:, 1), x(:, 3));
end

function [value, isterminal, direction] = eventfcn(t, x)
    boundary0 = [5 4];
    boundary1 = [25 20];
    value = [(x(1) - boundary0(1))*(boundary1(1) - x(1)), ...
             (x(2) - boundary0(2))*(boundary1(2) - x(2))];
    isterminal = [1, 1]; 
    direction = [0, 0];
end