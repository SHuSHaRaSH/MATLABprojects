function fres = f(t)
    fres = fillmissing(t.*cos(log(abs(t)))  , 'constant', 0);
end