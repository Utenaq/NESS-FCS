function [ corr2 ] = autocorrelation_ranged2( N)
% ������غ����ļ��㣬NΪtrace��corr_deltaΪ��غ���
[m,n] = size(N);
if(m==1)&&(n>1)
    N = N';
    m = n;
    n = 1;
end
p = zeros(224,1);
for i=1:32
    p(i) = i;
end
for j=1:12
    for i=1:16
        p(32+(j-1)*16+i) = p(32+(j-1)*16+i-1)+2^j;
    end
end
%t���ϵ�ȡ�㣬����ؿ��㷨һ��
mean = sum(N)/m;
dN = N-mean;
corr = zeros(224,1);
for i=1:224
    dN2 = zeros(m,1);
    dN2(1:m-p(i),1) = dN(p(i)+1:m,1);
    corr(i) = sum(dN2.*dN)/(m-p(i));
end
corr2 = corr/(mean^2);
end