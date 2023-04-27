%% functiion for Subtask 6
function dydx = linsyst(t, xy, a)

    if(size(a,1)==2 && size(a,2)==2)
        dydx = [(a(1, 1)*xy(1) + a(1, 2)*xy(2)); (a(2, 1)*xy(1) + a(2, 2)*xy(2))];
    else
        error('The equations degree is two');
    end
end
