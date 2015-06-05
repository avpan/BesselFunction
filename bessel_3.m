clear all;clc;
i = 16; %n terms
 j=1;
%finding error in bessel functions
for N = 1:2:i

    n(j) = N;
    z1 = N*pi/2;
    z2 = N*pi/4;
    for k = 0:20
        K(k+1) = k+1;
        I0(k+1,j) = (z1^2/4)^k/(factorial(k)*gamma(k+1));
        I1(k+1,j) = (z2^2/4)^k/(factorial(k)*gamma(k+2));
        Ineg1(k+1,j) = (z2^2/4)^k/(factorial(k)*gamma(k));
        
        besselI0(j) = sum(I0(:,j));
        besselI02(j) = integral(@(x) exp(z1*cos(x)),0,pi)/pi;
        besselI1(j) = sum(I1(:,j));
        besselIneg1(j) = sum(Ineg1(:,j));

       besselK02(j) = integral(@(x) cos(z1*x)./sqrt(x.^2+1),0,Inf);
        besselK0(j) = integral(@(x) exp(-z1*cosh(x)),0,Inf);
        besselK1(j) = gamma(1+1/2)*(2*z2)/sqrt(pi)*integral(@(x) ...
            cos(x)./(x.^2+z2.^2).^(3/2),0,Inf);
        
        %errors compared to matlab's built in bessel functions
        ErrI0(k+1,j) = abs((besseli(0,z1)-besselI0(j))/besseli(0,z1));
        ErrI1(k+1,j) = abs((besseli(1,z2)-besselI1(j))/besseli(1,z2));
        ErrK02(k+1,j) = abs((besselk(0,z1)-besselK02(j))/besselk(0,z1));
        ErrK0(k+1,j) = abs((besselk(0,z1)-besselK0(j))/besselk(0,z1));
        ErrK1(k+1,j) = abs((besselk(1,z2)-besselK1(j))/besselk(1,z2));
    end
    j = j+1;
end

% figure(1);
% subplot(2,2,1);
% plot(K,ErrI0(:,N),'black.');
% 
% subplot(2,2,2);
% plot(K,ErrI1(:,N),'black.');
% 
% subplot(2,2,3);
% plot(K,ErrK0(:,N),'black.');
% 
% subplot(2,2,4);
% plot(K,ErrK1(:,N),'black.');

[n,K] = meshgrid(n,K);
figure(1);
subplot(3,2,1);
plot3(K,n,ErrI0,'black.');
xlabel('k terms');ylabel('n terms');zlabel('error');
grid on;
title('Modifed Bessel Function: I_0(n\pi/2)');

subplot(3,2,2);
plot3(K,n,ErrI1,'black.');
xlabel('k terms');ylabel('n terms');zlabel('error');
grid on;
title('Modifed Bessel Function: I_1(n\pi/4)');

subplot(3,2,3);
plot3(K,n,ErrK02,'black.');
xlabel('k terms');ylabel('n terms');zlabel('error');
grid on;
title('Modifed Bessel Function: K0_{math}(n\pi/2)');

subplot(3,2,4);
plot3(K,n,ErrK0,'black.');
xlabel('k terms');ylabel('n terms');zlabel('error');
grid on;
title('Modifed Bessel Function: K0_{bell}(n\pi/2)');

subplot(3,2,5);
plot3(K,n,ErrK1,'black.');
xlabel('k terms');ylabel('n terms');zlabel('error');
grid on;
title('Modifed Bessel Function: K_1(n\pi/4)');
