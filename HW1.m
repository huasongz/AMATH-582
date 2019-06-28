clear all; close all; clc;
load Testdata
L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);


Uave = zeros(n,n,n);
for j = 1:20
    Utn(:,:,:) = fftn(reshape(Undata(j,:,:),n,n,n));
    Uave = Uave + Utn;
end
Uave = fftshift(Uave)./20;
[val,idx] = max(Uave(:));
[a,b,c] = ind2sub(size(Uave),idx)
isosurface(Kx,Ky,Kz,abs(Uave)./max(abs(Uave(:))),0.4)
xlabel('Kx')
ylabel('Ky')
zlabel('Kz')
title('frequency')
xx = Kx(a,b,c)
yy = Ky(a,b,c)
zz = Kz(a,b,c)
CF = sprintf('The central frequency in frequency domain is at %s %d %f.',xx,yy,zz);
disp(CF)


A = [];
B = [];
C = [];
filter=exp(-0.2*(((Kx-xx).^2)+((Ky-yy).^2)+((Kz-zz).^2)));
for i = 1:20
    dat(:,:,:) = reshape(Undata(i,:),n,n,n);
    Un = fftn(dat);
    Unt = fftshift(Un);
    Unft = filter.*Unt;
    Unf = ifftn(Unft);
    
    
    [val,idx] = max(Unf(:));
    [b,a,c] = ind2sub(size(Unf),idx);
    A(i) = a;
    B(i) = b;
    C(i) = c;
end
figure(2)
plot3(x(A),y(B),z(C))
xlabel('x')
ylabel('y')
zlabel('z')
title('Marble Movement Trajectory')
x20 = x(A(end));
y20 = y(B(end));
z20 = z(C(end));
L = sprintf('The 20th location of the marble is at %s %d %f.',x20,y20,z20);
disp(L)

