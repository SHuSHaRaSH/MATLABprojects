%% function for Subtask 10
function f = func1(t)
    f = @(t) cos(3*t).*(t>-1/2).*(t<1/2);
end