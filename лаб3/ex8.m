function ball
clear
clc

% начальные данные, ограничения и переменные ------------------------------
    t0 = 0;
    tf = 50;
    x0 = [7; 9];
    dx_dt0 = [-1; -1];
    y0 = [x0; dx_dt0];
    boundary0 = [5 4];
    boundary1 = [10 10];
    eps = 0.0000001;
    options = odeset('Events', @events);
    
% фигуры для графиков -----------------------------------------------------
    figure(1)
    cla
    hold on
    figure(2)
    cla
    hold on
    
% где будем хранить результаты численных решений (точки траектории и соотв время)
    tout = t0;
    yout = y0.';
    
% цикл вычисления точек траектории и времени, с учетом ограничений --------
    while(t0 < tf)
    %         disp('start');
        [t, y, te, ye, ie] = ode45(@f, [t0 tf], y0, options);

        tout = [tout; t];
        yout = [yout; y];
   
        t0 = te;
        if((t(end) >= tf) || (isempty(ie) == 1))
            break;
        end;
        % может выйти так, что за время вычисления траектории знаяения
        % несколько раз вышли за ограничения, поэтому делаем следующую
        % проверку
        if(size(ye, 1) == 1)
            % если вышли за или дошли до лев -ый/-ого кра -й/-я, то
            if((ye(1) - boundary0(1)) <= eps && ye(3) < 0)
                y0(3) = -ye(3);
                y0(1) = boundary0(1) + eps;
            % если вышли за или дошли до прав -ый/-ого кра -й/-я, то
            else if((boundary1(1) - ye(1)) <= eps && ye(3) > 0)
                y0(3) = -ye(3);
                y0(1) = boundary1(1) - eps;
            else
                y0(3) = ye(3);
                y0 (1) = ye(1);
                end;
            end;
            % если вышли за или дошли до нижн -ий/-его кра -й/-я, то
            if((ye(2) - boundary0(2)) <= eps && ye(4) < 0)
                y0(4) = -ye(4);
                y0(2) = boundary0(2) + eps;
            % если вышли за или дошли до верхн -ий/-его кра -й/-я, то
            else if((boundary1(2) - ye(2)) <= eps && ye(4) > 0)
                y0(4) = -ye(4);
                y0(2) = boundary1(2) - eps;
            else
                y0(4) = ye(4);
                y0(2) = ye(2);
                end;
            end;
        % несколько раз вышли за ограничения (все то же, только с учетом того, что ограничения сработали несколько раз)
        else
            if((ye(end, 1) - boundary0(1)) <= eps && ye(end, 3) < 0)
                y0(3) = -ye(end, 3);
                y0(1) = boundary0(1) + eps;
            else if((boundary1(1) - ye(end, 1)) <= eps && ye(end, 3) > 0)
                y0(3) = -ye(end, 3);
                y0(1) = boundary1(1) - eps;
            else
                y0(3) = ye(end, 3);
                y0(1) = ye(end, 1);
                end;
            end;
            if((ye(end, 2) - boundary0(2)) <= eps && ye(end, 4) < 0)
                y0(4) = -ye(end, 4);
                y0(2) = boundary0(2) + eps;
            else if((boundary1(2) - ye(end, 2)) <= eps && ye(end, 4) > 0)
                y0(4) = -ye(end, 4);
                y0(2) = boundary1(2) - eps;
            else
                y0(4) = ye(end, 4);
                y0(2) = ye(end, 2);
                end;
            end;
        end;
    end;
    
% график зависимости координат х и у от времени ---------------------------
    figure(1)
    plot(tout, yout(:, 1), 'r', tout, yout(:, 2), 'g');
    xlabel('time');
    ylabel('переменные х и у');
    legend('x = x(t)', 'y = y(t)');
    title('Зависимость переменных х и у от времени');
    hold off
    
% траектория мяча в промежутке времени [t0, tf] ---------------------------
    figure(2)
    %comet(yout(:, 1),yout(:, 2), 0.001);
    %plot(yout(:, 1), yout(:, 2), 'b');
    n = size(yout, 1);
    mov(1:n) = struct('cdata', [], 'colormap', []);
   for i = 1:n
       
       viscircles([yout(i,1), yout(i,2)], 0.5);
       ylim([3 20]);
       xlim([4 18]);
       
       daspect([1 1 1]);
       mov(i) = getframe();
       clf
   end
   movie(mov)
    
    xlabel('x');
    ylabel('y');
    title('Траектория мяча в плоскости, ограниченной перегородками');
    hold off
end

% система, которую решаем -------------------------------------------------
function dxdt = f(t, x)
    alpha = 0;
    dxdt = [x(3); x(4); -alpha*x(1); -alpha*x(2)];
end

% ограничения, которые передаются в решатель ОДУ ode45 --------------------
function [value, isterminal, direction] = events(t, x)
    boundary0 = [5 4];
    boundary1 = [25 20];
    value = [(x(1) - boundary0(1))*(boundary1(1) - x(1)), ...
             (x(2) - boundary0(2))*(boundary1(2) - x(2))];
    isterminal = [1, 1]; 
    direction = [0, 0];
end