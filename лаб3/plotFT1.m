%% function for Subtask 10
function S = plotFT(hFigure, fHandle, fFTHandle, step, inpLimVec, outLimVec)
    hold on;
    dist = inpLimVec(2)-inpLimVec(1);        
    N = round(dist/step);
    dt = dist/N;
    if(inpLimVec(1) > 0 || inpLimVec(2) < 0)
       inpLimVec = inpLimVec-dist*floor(inpLimVec(2)/dist);
    end
    t1 = 0:dt:inpLimVec(2)-dt;
    t2 = inpLimVec(1):dt:-dt;
    t = [t1,t2];
    t = t + 0.00001;
    
    f = fHandle(t);
    ft = fft(f);
    ft = ft*dt;
    FT = [ft,ft,ft,ft,ft,ft,ft,ft,ft,ft];
    
    if(~isempty(fFTHandle))
        w = linspace(outLimVec(1), outLimVec(2),100000);
        myft = fFTHandle(w);
    end
    df = 2*pi/dist;
    DF = -5*(N)*df:df:5*(N-1)*df+4*df;
    
    
    subplot(2,1,1);
    plot(DF,real(FT));
    if(~isempty(fFTHandle))
        hold on
        plot(w,real(myft));
    end
    axis([outLimVec(1), outLimVec(2), min(real(ft)), max(real(ft))]);
    legend({'Matlab''s', 'Mine'});
    subplot(2,1,2);
    plot(DF,imag(FT));
    if(~isempty(fFTHandle))
        hold on
        plot(w,imag(myft));
    end
    axis([outLimVec(1), outLimVec(2), min(imag(ft)), max(imag(ft))]);
    legend({'Matlab''s', 'Mine'});
    hold off
    
    S.nPoints = N;
    S.step = dt;
    S.inpLimVec = inpLimVec;
    S.outLimVec = outLimVec;    
end