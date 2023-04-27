0 = 0;
tf = 50;
x0 = [11; 11]; %starting point and upper boundary
dx_dt0 = [-1; -1];
y0 = [x0; dx_dt0];
boundary0 = [5 4];
boundary1 = [11 11];
eps = 0.0000001;

options = odeset('Events', @events);

alpha = 0;
f = @(t, x) [x(3); x(4); -alpha*x(1); -alpha*x(2)];
    
tout = t0; %time
yout = y0.'; %trajectory
    
while(t0 < tf) %getting points
    [t, y, t_curr, y_curr, ie] = ode45(f, [t0 tf], y0, options);
    tout = [tout; t];
    yout = [yout; y]; 
    t0 = t_curr;
    if((t(end) >= tf) || (isempty(ie) == 1))
        break;
    end
    if(size(y_curr, 1) == 1) %check that inside the boundaries
        if((y_curr(1) - boundary0(1)) <= eps && y_curr(3) < 0) %left
            y0(3) = -y_curr(3);
            y0(1) = boundary0(1) + eps;
        else
            if((boundary1(1) - y_curr(1)) <= eps && y_curr(3) > 0) %right
               y0(3) = -y_curr(3);
               y0(1) = boundary1(1) - eps;
            else
                y0(3) = y_curr(3);
                y0 (1) = y_curr(1);
            end
        end
        if((y_curr(2) - boundary0(2)) <= eps && y_curr(4) < 0) %lower
            y0(4) = -y_curr(4);
            y0(2) = boundary0(2) + eps;
        else
            if((boundary1(2) - y_curr(2)) <= eps && y_curr(4) > 0) %upper
                y0(4) = -y_curr(4);
                y0(2) = boundary1(2) - eps;
            else
                y0(4) = y_curr(4);
                y0(2) = y_curr(2);
            end
        end
    else %several boarders
            if((y_curr(end, 1) - boundary0(1)) <= eps && y_curr(end, 3) < 0)
                y0(3) = -y_curr(end, 3);
                y0(1) = boundary0(1) + eps;
            else
                if((boundary1(1) - y_curr(end, 1)) <= eps && y_curr(end, 3) > 0)
                    y0(3) = -y_curr(end, 3);
                    y0(1) = boundary1(1) - eps;
                else
                    y0(3) = y_curr(end, 3);
                    y0(1) = y_curr(end, 1);
                end
            end
            if((y_curr(end, 2) - boundary0(2)) <= eps && y_curr(end, 4) < 0)
                y0(4) = -y_curr(end, 4);
                y0(2) = boundary0(2) + eps;
            else
                if((boundary1(2) - y_curr(end, 2)) <= eps && y_curr(end, 4) > 0)
                    y0(4) = -y_curr(end, 4);
                    y0(2) = boundary1(2) - eps;
                else
                    y0(4) = y_curr(end, 4);
                    y0(2) = y_curr(end, 2);
                end
            end
        end
    end
    
    
    figure(1)
    plot(tout, yout(:, 1), 'r', tout, yout(:, 2), 'b');
    xlabel('t');
    ylabel('x / y');
    legend('x = x(t)', 'y = y(t)');
    title('Dependence of variables x y on time');
    
    figure(2)
    mov(1:length(yout)) = struct('cdata',[],'colormap',[]);
    hold on;
    rectangle('Position', [boundary0(1) boundary0(2) x0(1)-boundary0(1) x0(2)-boundary0(2)]);
    for i = 1: length(yout)
        plot(yout(1:i, 1), yout(1:i, 2), 'b');
        xlim([boundary0(1)-1  x0(1) + 1]);
        ylim([boundary0(2)-1 x0(2)+1]);
        mov(i) = getframe();
    end
    xlabel('x');
    ylabel('y');
    title('The trajectory of the ball');
    %movie(mov);
    
    clear;

function [value, isterminal, direction] = events(t, x)
    boundary0 = [5 4];
    boundary1 = [13 13];
    value = [(x(1) - boundary0(1))*(boundary1(1) - x(1)), ...
             (x(2) - boundary0(2))*(boundary1(2) - x(2))];
    isterminal = [1, 1]; 
    direction = [0, 0];
end