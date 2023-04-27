function dxydt = odefcn(t, xy)
    alpha = 1;
    dxydt = zeros(4,1);
    dxydt(1) =  3;
    dxydt(2) = -alpha* xy(1);
    dxydt(3) = -alpha*xy(4);
    dxydt(4) = 3;
end