function s = overlay(hFigure, fHandle, fFTHandle, step, inpLimVec, outLimVec)
    if ~isempty(hFigure)
        clf(hFigure)
    end
    
    info = get(hFigure, 'UserData');
    
    if ~isempty(info)
        if (isempty(outLimVec) == 1)
            outLimVec = get(info.real_ax, 'XLim');
        end
    else
        info = struct('real_axises', [], 'imag_axises', []);
    end
    
    g = @(t) fHandle(t + inpLimVec(1));
    T = inpLimVec(2) - inpLimVec(1);
    N = floor(T/step);
    x = linspace(0, T, N+1);
    x(end) = [];
    y = g(x);
    for i = 1:size(y, 2)
        if isnan(y(i))
            y(i) = y(i - 1);
        end
    end
    y = step*fft(y);
    sh = zeros(1, N);
    k0 = fix((-inpLimVec(1))/T*N); 
    for m = 1:N
        sh(m) = exp(-2*pi*1i*k0*(m-1)/N);
    end
    Y = y.*sh;
    %%step2 = 2*pi/T;
    step2 = 2*pi/T;
    L = step2*N;
    left = 0;
    right = L;
    nPer = 1;
    
    if isempty(info) & isempty(outLimVec)
        outLimVec = [0 L];
    end
    
    if(outLimVec(1) < left)
        c = ceil((left - outLimVec(1))/L);
        left = left - c*L;
        nPer = nPer + c
    end
    if(outLimVec(2) > right)
        c = ceil((outLimVec(2) - right)/L);
        right = right + c*L;
        nPer = nPer + c;
    end
    
    Y = repmat(Y, 1, nPer);
    X = left:step2:right;
    X(end) = [];
    NN = max((right - left)/step2 - 1, 20000);
    Xa = linspace(left, right, NN);
    
    
    pos1 = [0.1 0.15 0.35 0.7];
   %subplot('Position',pos1)
    if (~isempty(fFTHandle)) 
        FH = real(fFTHandle(Xa))
        realmin = min([real(Y), FH]);
        realmax = max([real(Y), FH]);
       
        legend("Численно(fft)", "Аналитически");
    else
        
         step3 = 20;
         Yo = zeros(1, N * nPer);
         hold on
        for i = -1:1
            Yi = circshift(Y,i*step3)
            p1 = plot(X, real(Yi),'r');
            
            p1.LineWidth = 1;
            disp(size(Yi));
            disp(size(Yo));
            Yo = Yi + Yo; 
            
            
        end
        
        realmin = min(real(Yo));
        realmax = max(real(Yo));
        p2 = plot(X, real(Yo),'gr');
        p2.LineWidth = 1;
        legend([p1 p2],{"исходный","наложение"});
    end
    ylabel('real part of F');
    txt = texlabel('lambda'); 
    xlabel(txt)
    if abs(realmin) + abs(realmax) 
        axis([outLimVec(1) outLimVec(2) realmin realmax]);
    else
        axis([outLimVec(1) outLimVec(2) -0.01 0.01]);
    end
    txt = texlabel('Real F'); 
    title(txt)
    info.real_axises = gca;
    exportgraphics(hFigure,'spectrumOverlay3.eps')
    hold on
    pos2 = [0.6 0.15 0.35 0.7];
  % subplot('Position',pos2)
    imY = imag(Y);
%     if (~isempty(fFTHandle))
%         FH = imag(fFTHandle(Xa));
%         imagmin = min([imY, FH]);
%         imagmax = max([imY, FH]);
%         
%         %%p1 = plot(X, imY);
%         p1.LineWidth = 1;
%         hold on
%         %%p2 = plot(Xa, FH);
%         p2.LineWidth = 1;
%         legend("Численно(fft)", "Аналитически");
%     else
%         imagmin = min(imY);
%         imagmax = max(imY);
%         %%plot(X, imY);
%         legend("Численно(fft)");
%     end
    ylabel('imag part of F');
    txt = texlabel('lambda'); 
    xlabel(txt)
    if abs(imagmin) + abs(imagmax) 
        axis([outLimVec(1) outLimVec(2) imagmin imagmax]);
    else
        axis([outLimVec(1) outLimVec(2) -0.01 0.01]);
    end
    title('Imag F');
    
    info.imag_axises = gca;
    
    set(hFigure, 'UserData', info);
    s = struct('nPoints', N, 'step', step, 'inpLimVec', inpLimVec, 'outLimVec', outLimVec);
    
    
end