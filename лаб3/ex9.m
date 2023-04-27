function ex9()
%ЗАДАНИЕ:
% Рассмотреть систему двух тел на плоскости:
% m1*x1 = G*m1*m2*(x2 ? x1) / norm(x1-x2)^3, x1 ? R2
% m2*x2 = G*m1*m2*(x1 ? x2) / norm(x1-x2)^3, x2 ? R2
% Решить систему численно . Нарисовать на плоскости анимацию движения 
% траекторий x1(t), x2(t). Подобрать параметры так, чтобы 
% продемонстрировать движение двух типов: по пересекающимся орбитам 
% («восьмёрка») и вокруг общего центра.

% очистка -----------------------------------------------------------------
clear
clc

% начальные данные, переменные и ограничения ------------------------------
    % эти начальные данные красивее, чем, например, при [0; 0] [5; 5] [0; 25] [-6; 0]
    x10 = [0; 0];
    x20 = [0; 25];
    dx1_dt0 = [5; 5];
    dx2_dt0 = [10; 0];
    y0 = [x10; x20; dx1_dt0; dx2_dt0];
    
    t0 = 0;
    t1 = 15;
    
% решаем систему и рисуем решение -----------------------------------------
    % ---------------------- example 1 ----------------------------
    answer = input('Если хотите посмотреть пример при G = 6.67 * 10^(-2), m1 = 400000, m2 = 100 нажмите 1, иначе нажмите 0');
    if(answer == 1)
        [t, y] = ode45(@f1, [t0, t1], y0);
        % график и анимация
        figure(1)
        cla
        hold on
        title('G = 6.67 * 10^(-3), m1 = 100000000, m2 = 100');
        n = min(size(t, 1), 3000);
        for i = 1:n
            plot(y(1:i, 1), y(1:i, 2),'c', y(1:i, 3), y(1:i, 4), 'm');
            pause(0.01);
        end;
        hold off
    end;
    
    % ---------------------- example 2 ----------------------------
    answer = input('Если хотите посмотреть пример при G = 6.67 * 10^(-2), m1 = 400000, m2 = 100 нажмите 1, иначе нажмите 0');
    if(answer == 1)
        [t, y] = ode45(@f2, [t0, t1], y0);
        % график и анимация
        figure(2)
        cla
        hold on
        title('G = 6.67 * 10^(-2), m1 = 400000, m2 = 100');
        n = min(size(t, 1), 1200);
        for i = 1:n
            plot(y(1:i, 1), y(1:i, 2),'c', y(1:i, 3), y(1:i, 4), 'm');
            pause(0.01);
        end;
        hold off
    end;
    
    % ---------------------- example 3 ----------------------------
    answer = input('Если хотите посмотреть пример при G = 6.67, m1 = 10000000, m2 = 100 нажмите 1, иначе нажмите 0');
    if(answer == 1)
        [t, y] = ode45(@f3, [t0, t1], y0);
        % график и анимация
        figure(3)
        cla
        hold on
        title('G = 6.67, m1 = 10000000, m2 = 100');
        n = min(size(t, 1), 1500);
        for i = 1:n
            plot(y(1:i, 1), y(1:i, 2),'c', y(1:i, 3), y(1:i, 4), 'm');
            pause(0.01);
        end;
        hold off
    end;
    
    % ---------------------- example 4 ----------------------------
    answer = input('Если хотите посмотреть пример при G = 3.67 * 10^(-1), m1 = 120, m2 = 100 нажмите 1, иначе нажмите 0');
    if(answer == 1)
        % эти начальные данные красивее, чем, например, при [0; 0] [5; 5] [0; 25] [-6; 0]
        x10 = [0; 0];
        x20 = [0; 0.4];
        dx1_dt0 = [5; 5];
        dx2_dt0 = [-5; -5];
        y0 = [x10; x20; dx1_dt0; dx2_dt0];

        t0 = 0;
        t1 = 10;

        [t, y] = ode45(@f4, [t0, t1], y0);
        % график и анимация
        figure(2)
        cla
        hold on
        title('G = 3.67 * 10^(-1), m1 = 120, m2 = 100');
        n = min(size(t, 1), 1200);
        for i = 1:n
            plot(y(1:i, 1), y(1:i, 2),'c', y(1:i, 3), y(1:i, 4), 'm');
            pause(0.01);
        end;
        hold off
    end;
    % ---------------------- example 5 ----------------------------
    answer = input('Если хотите посмотреть пример при G = 6.67 * 10^(-1), m1 = 220, m2 = 100 нажмите 1, иначе нажмите 0');
    if(answer == 1)
        % эти начальные данные красивее, чем, например, при [0; 0] [5; 5] [0; 25] [-6; 0]
        x10 = [0; 0];
        x20 = [0; 0.4];
        dx1_dt0 = [5; 5];
        dx2_dt0 = [-5; -5];
        y0 = [x10; x20; dx1_dt0; dx2_dt0];

        t0 = 0;
        t1 = 10;

        [t, y] = ode45(@f5, [t0, t1], y0);
        % график и анимация
        figure(2)
        cla
        hold on
        title('G = 6.67 * 10^(-1), m1 = 20, m2 = 100');
        n = min(size(t, 1), 1200);
        for i = 1:n
            plot(y(1:i, 1), y(1:i, 2),'c', y(1:i, 3), y(1:i, 4), 'm');
            pause(0.01);
        end;
        hold off
    end;
    
end

% система для первого примера ---------------------------------------------
function dxdt = f1(t, x)
    
    G = 6.67 * 10^(-3);
    m1 = 100000000;
    m2 = 100;
    
    x11 = x(1);
    x12 = x(2);
    x21 = x(3);
    x22 = x(4);
    dx11_dt = x(5);
    dx12_dt = x(6);
    dx21_dt = x(7);
    dx22_dt = x(8);
    
    norm = @(x) sqrt(x(1)^2 + x(2)^2);
    n12_3 = norm([x11 - x21; x12 - x22])^3;
    
    d2x11_dt2 = G * m2 * (x21 - x11) / n12_3;
    d2x12_dt2 = G * m2 * (x22 - x12) / n12_3;
    d2x21_dt2 = G * m1 * (x11 - x21) / n12_3;
    d2x22_dt2 = G * m1 * (x12 - x22) / n12_3;
    
    dxdt = [dx11_dt;    dx12_dt;    dx21_dt;    dx22_dt; ...
            d2x11_dt2;  d2x12_dt2;  d2x21_dt2;   d2x22_dt2];
end

% система для второго примера ---------------------------------------------
function dxdt = f2(t, x)
    
    G = 6.67 * 10^(-2);
    m1 = 400000;
    m2 = 100;
    
    x11 = x(1);
    x12 = x(2);
    x21 = x(3);
    x22 = x(4);
    dx11_dt = x(5);
    dx12_dt = x(6);
    dx21_dt = x(7);
    dx22_dt = x(8);
    
    norm = @(x) sqrt(x(1)^2 + x(2)^2);
    n12_3 = norm([x11 - x21; x12 - x22])^3;
    
    d2x11_dt2 = G * m2 * (x21 - x11) / n12_3;
    d2x12_dt2 = G * m2 * (x22 - x12) / n12_3;
    d2x21_dt2 = G * m1 * (x11 - x21) / n12_3;
    d2x22_dt2 = G * m1 * (x12 - x22) / n12_3;
    
    dxdt = [dx11_dt;    dx12_dt;    dx21_dt;    dx22_dt; ...
            d2x11_dt2;  d2x12_dt2;  d2x21_dt2;   d2x22_dt2];
end

% система для третьего примера --------------------------------------------
function dxdt = f3(t, x)
    
    G = 6.67;
    m1 = 1000000;
    m2 = 100;
    
    x11 = x(1);
    x12 = x(2);
    x21 = x(3);
    x22 = x(4);
    dx11_dt = x(5);
    dx12_dt = x(6);
    dx21_dt = x(7);
    dx22_dt = x(8);
    
    norm = @(x) sqrt(x(1)^2 + x(2)^2);
    n12_3 = norm([x11 - x21; x12 - x22])^3;
    
    d2x11_dt2 = G * m2 * (x21 - x11) / n12_3;
    d2x12_dt2 = G * m2 * (x22 - x12) / n12_3;
    d2x21_dt2 = G * m1 * (x11 - x21) / n12_3;
    d2x22_dt2 = G * m1 * (x12 - x22) / n12_3;
    
    dxdt = [dx11_dt;    dx12_dt;    dx21_dt;    dx22_dt; ...
            d2x11_dt2;  d2x12_dt2;  d2x21_dt2;   d2x22_dt2];
end

% система для четвертого примера --------------------------------------------
function dxdt = f4(t, x)
    
    G = 3.67 * 10^(-1);
    m1 = 120;
    m2 = 100;
    
    x11 = x(1);
    x12 = x(2);
    x21 = x(3);
    x22 = x(4);
    dx11_dt = x(5);
    dx12_dt = x(6);
    dx21_dt = x(7);
    dx22_dt = x(8);
    
    norm = @(x) sqrt(x(1)^2 + x(2)^2);
    n12_3 = norm([x11 - x21; x12 - x22])^3;
    
    d2x11_dt2 = G * m2 * (x21 - x11) / n12_3;
    d2x12_dt2 = G * m2 * (x22 - x12) / n12_3;
    d2x21_dt2 = G * m1 * (x11 - x21) / n12_3;
    d2x22_dt2 = G * m1 * (x12 - x22) / n12_3;
    
    dxdt = [dx11_dt;    dx12_dt;    dx21_dt;    dx22_dt; ...
            d2x11_dt2;  d2x12_dt2;  d2x21_dt2;   d2x22_dt2];
end

% система для четвертого примера --------------------------------------------
function dxdt = f5(t, x)
    
    G = 6.67 * 10^(-1);
    m1 = 200;
    m2 = 100;
    
    x11 = x(1);
    x12 = x(2);
    x21 = x(3);
    x22 = x(4);
    dx11_dt = x(5);
    dx12_dt = x(6);
    dx21_dt = x(7);
    dx22_dt = x(8);
    
    norm = @(x) sqrt(x(1)^2 + x(2)^2);
    n12_3 = norm([x11 - x21; x12 - x22])^3;
    
    d2x11_dt2 = G * m2 * (x21 - x11) / n12_3;
    d2x12_dt2 = G * m2 * (x22 - x12) / n12_3;
    d2x21_dt2 = G * m1 * (x11 - x21) / n12_3;
    d2x22_dt2 = G * m1 * (x12 - x22) / n12_3;
    
    dxdt = [dx11_dt;    dx12_dt;    dx21_dt;    dx22_dt; ...
            d2x11_dt2;  d2x12_dt2;  d2x21_dt2;   d2x22_dt2];
end