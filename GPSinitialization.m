function GD = GPSinitialization(N, dim, lb, ub)

%% ����ѵ㼯
tmp1 = (1:N)'*ones(1, dim);
Ind = 1:dim;
prime1 = primes(100*dim);           % �ҳ�����
q = find(prime1 >= (2*dim+3));
tmp2 = (2*pi.*Ind)/prime1(q(1));
tmp2 = 2*cos(tmp2);
tmp2 = ones(N, 1)*tmp2;
GD = tmp1.*tmp2;
GD = mod(GD, 1);

GD = GD.*(ub-lb)+lb;