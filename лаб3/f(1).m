%% function for Subtask 4
function dxydt = f(t, xy)
    dxydt = zeros(4,1);
    dxydt(1) = xy(3);
    dxydt(2) = xy(4);
    dxydt(3) = 0;
    dxydt(4) = 0;
end