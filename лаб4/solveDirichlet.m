% function for Subtask 6
function unum = solveDirichlet(fHandle,xiHandle,etaHandle,mu,N,M)
    deltax = 1/% function for Subtask 6
function unum = solveDirichlet(fHandle,xiHandle,etaHandle,mu,N,M)
    deltax = 1/N;
    deltay = 1/M;
    x = linspace(0,1-deltax,N);
    y = linspace(0,1-deltay,M);
    [X,Y] = meshgrid(y,x);

    f = fHandle(X,Y);

    f(1:N,1)=0;
    f(1,1:M)=0;

    xi = zeros(1,M);
    xi(1:M) = xiHandle(y(1:M));

    eta = zeros(1,N);
    eta(1:N) = etaHandle(x(1:N));

    alpha = ifft(eta);
    beta = ifft(xi);
    C = zeros(N,M);
    for p=1:N
        C(p,1:M) = -4*(sin(pi*(p-1)/N)).^2/(deltax.^2)-4*(sin(pi*((1:M)-1)/M)).^2/(deltay.^2)-mu;
    end
    Q = zeros(N+M-1,N+M-1);
    t = zeros(1,N+M-1);
    %��������� �������� N �������
    Dmat = ifft2(f);
    for p = 1:N
        % �������� �������� ������ �����
        D = sum((1./C(p,1:N)).*Dmat(p,1:N));
        t(p) = alpha(p)-D;

        Cdivided=1./C(p,1:M);

        Yl = ifft(Cdivided)/N;
        % ������ ������� ��������� ��������
        Q(p,2:M)=Yl(2:M);
        % ��������� ������ M ��������� p-��� �������

        Q(p,1) = Q(p,1)+sum(Cdivided)./(M*N);
        % ������ ������� ��������� ��������


        Q(p,(2:N)+M-1) = (sum(Cdivided)./(M*N))*exp((2*pi*1j*((2:N)-1)*(p-1))/N);
        %��������� M+1 ...N+M-1 ����� �������

        % ��������� (p,1) ��������� ���������
    end

    %��������� ������ ����� �������
    for q = 2:M
        D = sum((1./C(p,1:N)).*Dmat(p,1:N));
        t(q+N-1) = beta(q)-D;
        % ��������� ������ �����

        Cdivided(1:N)=1./C(1:N,q);

        %�������� ��������������� ��������

        Q(q+N-1,1) = Q(q+N-1,1)+sum(Cdivided)./(M*N);
        % ��������� � 0(1) �������

        Q(q+N-1,2:M) = (sum(Cdivided)./(M*N))*exp((2*pi*1j*((2:M)-1)*(q-1))/M);
        % ��������� ������ M ���������

        Xk = ifft(Cdivided)/M;

        %��������� � 0(1) ��������
        Q(q+N-1,M+(2:N)-1) = Xk(2:N);
    end
    f_NM = bicg(Q,t',1e-7,1000);
    % ������ �������
    f(1,1:M) = f_NM(1:M);
    f(2:N,1) = f_NM(M+(2:N)-1);
    % ��������� ������� f ������������ ����������
    bpq = ifft2(f);
    apq = bpq./C;
    unum = real(fft2(apq));
end