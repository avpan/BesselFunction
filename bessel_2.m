% clear all;clc;
% i = 15; %n terms
project;

%This is to construct the modified bessel functions manually in order to
%compare with MATLAB's built in functions
I0 = 0;
I1 = 0;
Ineg1 = 0;
j = 1; 
for N = 1:2:i
    z1 = N*pi/2;
    z2 = N*pi/4;
    for k = 0:20
        I0(k+1,j) = (z1^2/4)^k/(factorial(k)*gamma(k+1));
        I1(k+1,j) = (z2^2/4)^k/(factorial(k)*gamma(k+2));
    end
    besselI0(j) = sum(I0(:,j));
    %besselI02(j) = integral(@(x) exp(z1*cos(x)),0,pi)/pi;
    besselI1(j) = sum(I1(:,j))*(z2/2);
    %besselI12(j) = integral(@(x) exp(z1*cos(x)).*cos(N*x),0,pi)/pi;
    besselK02(j) = integral(@(x) cos(z1*x)./sqrt(x.^2+1),0,Inf);
    besselK0(j) = integral(@(x) exp(-z1*cosh(x)),0,Inf);
    besselK1(j) = gamma(1+1/2)*(2*z2)/sqrt(pi).*integral(@(x)cos(x)./(x.^2+z2.^2).^(3/2),0,Inf);
    j = j+1;
 end


j=1;
for N = 1:2:i
    n(j) = N;
    k = N*pi/2;
    P = sin(k)/besselI0(j)*(besselI0(j)/(N*pi)-.25*(besselI0(j)...
        *besselK1(j) + besselK0(j)*besselI1(j)));
    if j == 1
        phi3(j) = P;
    else
        phi3(j) = P + phi3(j-1);
    end
    j = j+1;
end
phi3 = 2.*phi3;

j=1;
for N = 1:2:i
    n(j) = N;
    k = N*pi/2;
    P = sin(k)/besselI0(j)*(besselI0(j)/(N*pi)-.25*(besselI0(j)...
        *besselK1(j) + besselK02(j)*besselI1(j)));
    if j == 1
        phi4(j) = P;
    else
        phi4(j) = P + phi4(j-1);
    end
    j = j+1;
end
phi4 = 2.*phi4;

err = abs((phi-phi4)./phi);
err2 = abs((phi-phi3)./phi);

figure(1);
subplot(1,2,1);
plot(n,phi4,'black.');
xlabel('n terms');ylabel('\Phi/\Phi_{disk}');
title('Ratio using K0_{math}');
grid on;

figure(1);
subplot(1,2,2);
plot(n,phi3,'black.');
xlabel('n terms');ylabel('\Phi/\Phi_{disk}');
title('Ratio using K0_{bell}');
grid on;

figure(3);
% subplot(2,2,1);
% plot(n,err,'black.');
% xlabel('n terms');ylabel('error diff.');
% title('error diff. using K0_{math} (A)');
% grid on;

subplot(1,2,1);
plot(n(3:i/2),err(3:i/2),'black.');
xlabel('n terms');ylabel('error diff.');
title('error diff. using K0_{math} (B)');
grid on;
% 
% subplot(2,2,3);
% plot(n,err2,'black.');
% xlabel('n terms');ylabel('error diff.');
% title('error diff. sing K0_{bell} (A)');
% grid on;

subplot(1,2,2);
plot(n(3:i/2),err2(3:i/2),'black.');
xlabel('n terms');ylabel('error diff.');
title('error diff. using K0_{bell}(B)');
grid on;

% figure(4);
% subplot(2,2,1);
% plot(n,besselI0,'black.');
% xlabel('n terms');ylabel('I_0(n\pi/2)');
% title('M. Bessel Function I_0')
% grid on;
% 
% subplot(2,2,2);
% plot(n,besselI1,'black.');
% xlabel('n terms');ylabel('I_1(n\pi/4)');
% title('M. Bessel Function I_1')
% grid on;
% 
% subplot(2,2,3);
% plot(n,besselK0,'black.');
% xlabel('n terms');ylabel('K_0(n\pi/2)');
% title('M. Bessel Function K0_{bell}')
% grid on;

% subplot(2,2,3);
% plot(n,besselK1,'black.');
% xlabel('n terms');ylabel('K_1(n\pi/4)');
% title('M. Bessel Function K_1')
% grid on;
% 
% figure(5);
% subplot(1,2,1);
% plot(n,besselK02,'black.',n,besselMATK0,'blue*');
% xlabel('n terms');ylabel('K_0(n\pi/2)');
% title('M. Bessel Function K0_{math}')
% grid on;
% 
% subplot(1,2,2);
% plot(n,besselK0,'black.');
% xlabel('n terms');ylabel('K_0(n\pi/2)');
% title('M. Bessel Function K0_{bell}')
% grid on;

% figure(5);
% subplot(1,2,1);
% plot(n,besselI1.*besselK0,'black.');
% xlabel('n terms');ylabel('I_1(n\pi/4)K_0(n\pi/2)');
% title('M. Bessel Function K_0')
% grid on;
% 
% subplot(1,2,2);
% plot(n,besselK1.*besselI0,'black.');
% xlabel('n terms');ylabel('I_0(n\pi/2)K_1(n\pi/4)');
% title('M. Bessel Function K_1')
% grid on;

