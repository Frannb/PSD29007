% Filtro chebyshev 1 - HP (Passa-alta)
close all;
clear all;
clc;

% Dados do filtro
fs = 40;        % é a frequência de rejeição em Hz (stopband)
fp = 200;       % é a frequência de passagem em Hz
Ap = 1;         % é a atenuação em dB na frequência de passagem
As = 50;        % é a atenuação em dB na frequência de stopband
Gtopo = 10;

Ws = 2*pi*fs;
Wp = 2*pi*fp; 

% Transformações de frequência de filtros analógicos
Os = Wp/Ws;
Op = 1;

E = sqrt(10.^(Ap/10)-1); % Épsilon
% Determinando a ordem
n = ceil((acosh (sqrt((10.^(As/10)-1)/E^2)))/(acosh(Os)));

% Obtendo os polos do filtro:
k = 1:n; % k varia de 1 até n 
Tk = ((2*k-1)*pi)/(2*n);
F2 = (1/n)*asinh(1/E);
p = -sinh(F2)*sin(Tk)+1j*cosh(F2)*cos(Tk);

% Obtendo a função de transferência:
D = real(poly(p)); % Pega o denominador
H0 = D(end)/sqrt(1+(E^2)); % Ajustando H0
syms p s;
Hp(p) = H0 / (D(1)*p^4 + D(2)*p^3 + D(3)*p^2 + D(4)*p + D(5));
pretty(vpa(Hp,3))

% Plotando os plos 
figure(1);
zplane(1, D);
%plot(real(p), imag(p), 'x');
title('Polos do prototipo');
axis([-1.25 1.25 -1.25 1.25]);
grid on;

% Resposta em frequência do filtro passa baixa 
figure(2)
freqs(H0, D)
title('Resposta em frequência do filtro passa baixa');

Hs(s)= subs(Hp, Wp/s);
pretty(vpa(simplify(Hs),4))

% Obtendo numerador e denominador
[Ns, Ds] = numden(Hs);

Nsv = sym2poly(Ns);
Dsv = sym2poly(Ds);

% Normalizar
Nsvn = Nsv/Dsv(1) * 10^((Gtopo)/20);
Dsvn = Dsv/Dsv(1);

% Resposta em frequência do filtro passa alta 
figure(3)
[h, w] = freqs(Nsvn, Dsvn, logspace(2, 5, 10000));
semilogx(w,20*log10(abs(h)));
grid on;
hold on;

plot(Ws,Gtopo-As,'or')
plot(Wp,Gtopo-Ap,'or')
ylim([-50 12]);

line([1256.64 1256.64 98*10^3],[-50 9 9],'Color','red','LineStyle','--');
line([10^2 250 250],[-40 -40 12],'Color','red','LineStyle','--');
xlabel('x');
ylabel('y');

title(sprintf('Filtro passa-alta chebyshev tipo 1, n = %d',n))

figure(4);
zplane(1, Dsv);

%%
fs = 40;        % é a frequência de rejeição em Hz (stopband)
fp = 200;       % é a frequência de passagem em Hz
Ap = 1;         % é a atenuação em dB na frequência de passagem
As = 50;        % é a atenuação em dB na frequência de stopband
Gtopo = 10;

[n,Wn] = cheb1ord(Wp,Ws,Ap,As,'s') % n: ordem do filtro, wn: freq de corte do filtro
[b,a] = cheby1(n,Ap,Wp,'high', 's'); % função de transferencia 
b = b*10^(Gtopo/20);

figure(1);
[h,w] = freqs(b,a,logspace(2, 5, 10000));
semilogx(w,20*log10(abs(h)));
ylim([8 10.5]);
grid on;
title(sprintf('Filtro passa-alta chebyshev tipo 1, n = %d',n))


