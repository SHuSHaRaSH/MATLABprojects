clear;
 
time0 = 0;
time1 = 0.1;
 
init_cond = [40; 0; 0; 5 * -10^5; 50; 0; 0; 5 * 10^5];
[t,x] = ode45(@twoBody, [time0 time1], init_cond);
 
%axis ([min(min(x(:,1),x(:,5))) max(max(x(:,1),x(:,5))) min(min(x(:,2),x(:,6))) max(max(x(:,2),x(:,6)))])
axis([min(min(x(:,1),x(:,5))) max(max(x(:,1),x(:,5))) -30 30]);
hold on;
axis manual;
 
for i = 2:length(t)
    plot(x(i-1:i, 1), x(i-1:i, 2), ':ro');
    hold on;
    plot(x(i-1:i, 5), x(i-1:i, 6), ':bo');
    legend("body 1", "body 2");
    pause(0.001);
end
 
function dxdt = twoBody(t, x)
    G  = 6.67*10^(-11);  
    m1 = 11^22;
    m2 = 11^22;%%%%%%%%%%%%%%%%%%%%%%%%%%
    dxdt = [x(3); 
            x(4);
            G*m2*(x(5)-x(1))/sqrt((x(5)-x(1))^2+(x(6)-x(2))^2)^3;
            G*m2*(x(6)-x(2))/sqrt((x(5)-x(1))^2+(x(6)-x(2))^2)^3; 
            x(7); 
            x(8);
            -G*m1*(x(5)-x(1))/sqrt((x(5)-x(1))^2+(x(6)-x(2))^2)^3;
            -G*m1*(x(6)-x(2))/sqrt((x(5)-x(1))^2+(x(6)-x(2))^2)^3;
            ];
end
 