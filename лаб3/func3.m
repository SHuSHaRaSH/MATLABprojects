%% function for Subtask 10
function f = func3(t)
    f = @(t) exp(-abs(t).^3).*t.^5;
end