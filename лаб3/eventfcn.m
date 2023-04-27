function [value, isterminal, direction] = eventsfcn(t, x)
    boundary0 = [-5 5];
    boundary1 = [-10 10];
    value = [(x(1) - boundary0(1))*(boundary0(2) - x(1)), ...
             (x(3) - boundary1(1))*(boundary1(2) - x(3))];
    isterminal = [1, 1]; 
    direction = [0, 0];
end