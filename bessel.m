clear all;clc;
i = 500; %# of n terms

j = 1;
%double precision 64-bit
for N = 1:2:i %sum over odd terms
    n(j) = N;
    k = N*pi/2;
    P = sin(k)/besseli(0,k)*(besseli(0,k)/(N*pi)-.25*(besseli(0,k)*besselk(1,k/2)+...
        besselk(0,k)*besseli(1,k/2)));
    P = double(P); %matlabs default is double precision, but in case
    if j == 1
        phi(j) = double(P);
    else
        phi(j) = double(P + phi(j-1));
    end
    j = j+1;
end
phi = 2.*phi;

j=1;
%single precision 32-bit
for N = 1:2:i %sum over odd terms
    n(j) = N;
    k = N*pi/2;
    P = sin(k)/besseli(0,k)*(besseli(0,k)/(N*pi)-(1/4)*(besseli(0,k)*...
        besselk(1,k/2)+besselk(0,k)*besseli(1,k/2)));
    P = single(P); %convert to single precision
    if j == 1
        phi2(j) = single(P);
    else
        phi2(j) = single(P + phi2(j-1));
    end
    j = j+1;
end
phi2 = single(2.*phi2);
% figure(1);
% subplot(1,2,1);
% plot(n,phi,'black.');
% xlabel('n terms');ylabel('\Phi/\Phi_{disk}');
% title('double precision');
% grid on;
% subplot(1,2,2);
% plot(n,phi2,'black.');
% xlabel('n terms');ylabel('\Phi/\Phi_{disk}');
% title('single precision');
% grid on;
% 
% figure(2);
% plot(n,abs(phi-phi2),'black.');
% xlabel('n terms');ylabel('|\Phi_{dp}-\Phi_{sp}|/\Phi_{disk}');
% title('precision difference');
% grid on;
