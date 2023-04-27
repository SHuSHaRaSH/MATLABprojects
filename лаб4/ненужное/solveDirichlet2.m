function unum = solveDirichlet2(fHandle,xiHandle,etaHandle,mu,N,M)

    deltax = 1/N;
    deltay = 1/M;
    x = linspace(0,1function unum = solveDirichlet2(fHandle,xiHandle,etaHandle,mu,N,M)

    deltax = 1/N;
    deltay = 1/M;
    x = linspace(0,1-deltax,N);
    y = linspace(0,1-deltay,M);
    [X,Y] = meshgrid(y,x);

    f = fHandle(X,Y);
    
    
    % для поиска бкс
    f(1:N,1)=0;
    f(1,1:M)=0;

    xi(1:M) = xiHandle(y(1:M));

    eta(1:N) = etaHandle(x(1:N));

    gamma = ifft(xi);
    delta = ifft(eta);
    
    C = zeros(N,M);
    % находим С_кс
    for p=1:N
        C(p,1:M) = -4*(sin(pi*(p-1)/N)).^2/(deltax.^2)-4*(sin(pi*((1:M)-1)/M)).^2/(deltay.^2)-mu;
    end
    
    mainmatrix = zeros(N+M-1,N+M-1);
    mainvector = zeros(1,N+M-1);
    %заполняем перввыми N числами
    
    % B_ks
    Bks = ifft2(f);
    for p=1:N
        % А_к
        Ak = (1./C(p,1:M)).*Bks(p,1:M);
        % гамма
        D = sum(Ak);
        
        mainvector(p) = delta(p)-D;

        Cdivided(1:M)=1./C(p,1:M);
        %2
        tao = ifft(Cdivided)/N;
        % первый элемент заполняем отдельно
        mainmatrix(p,2:M)=tao(2:M);
        % Заполняем первые M элементов p-ого столбца

        mainmatrix(p,1) = mainmatrix(p,1)+sum(Cdivided)./(M*N);
        % первый элемент дополняем отдельно

        mainmatrix(p,(2:N)+M-1) = (sum(Cdivided)./(M*N))*exp((2*pi*j*((2:N)-1)*(p-1))/N);
        %Вычисляем M+1 ...N+M-1 коэфф системы

        % Добавляем (p,1) последнее слагаемое
    end

    %Заполняем вторую часть матрицы
    for q = 2:M
        Ak = (1./C(1:N,q)).*Bks(1:N,q);
        D = sum(Ak);
        mainvector(q+N-1) = gamma(q)-D;
        % вычисляем правую часть

        Cdivided(1:N)=1./C(1:N,q);

        %вычислям вспомогательные значения

        mainmatrix(q+N-1,1) = mainmatrix(q+N-1,1)+sum(Cdivided)./(M*N);
        % Добавляем к 0(1) элменту

        mainmatrix(q+N-1,2:M) = (sum(Cdivided)./(M*N))*exp((2*pi*1j*((2:M)-1)*(q-1))/M);
        % Заполняем первые M элементов

        dzeta = ifft(Cdivided)/M;

        
        mainmatrix(q+N-1,M+(2:N)-1) = dzeta(2:N);
    end
    
    fNM = cgs(mainmatrix,mainvector',1e-7,1000);
    % Решаем систему
    f(1,1:M) = fNM(1:M);
    f(2:N,1) = fNM(M+(2:N)-1);
    % Заполняем матрицу f полученнымим значениями
    bpq = ifft2(f);

    apq = zeros (N,M);
    apq = bpq./C;
    unum = real(fft2 (apq));
end
