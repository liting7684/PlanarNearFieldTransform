function [Amp]=Current_Excitation(M,N,dx,dy)%��������ֲ�Amn

midM=(M-1)/2;midN=(N-1)/2;
Xm=[];Yn=[];%��(m, n)�����е�Ԫ�����ĵ������
for m=1:1:M
    Xm(m)=(m-midM-1)*dx;
end
% ������yn(n)����������y������꣬��֤��˼·�ǵ�16������Ϊ0
for n=1:1:N
    Yn(n)=(n-midN-1)*dy;
end
ampx=ones(1,M);%��x����Ϊ�ȷ���
ampy=ones(1,N);
%y�������ҷֲ�
for i=1:N
    ampy(i)=cos(pi*(i-midN-1)/N);
end

% x�����б�ѩ��ֲ�
SLL=-55;%�����ƽ
ampx=Chebyshev_Dist_Odd(M,SLL);
% ampx=chebwin1(43,-55);

Amp=ampx'*ampy;
[ymesh,xmesh]=meshgrid(Yn,Xm);
figure(1);
surf(xmesh,ymesh,Amp);
xlabel('xdimension');ylabel('ydimension');title('�Գ����ӷ��ȷֲ�');
disp('__________Complete the calculation of current__________');
% disp(Amp);
end