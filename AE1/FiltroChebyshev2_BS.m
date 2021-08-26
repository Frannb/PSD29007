% https://wiki.sj.ifsc.edu.br/index.php/Fazer_mascara_de_um_filtro_no_Matlab
% Filtro chebyshev 2 - BS (rejeita-faixa) 
close all;
clear all;
clc;

% Dados do filtro
fp1 = 1200;
fs1 = 2000;
fs2 = 3000;
fp2 = 3800;
fa = 8000;
Ap = 3;           
As = 40;       
Gtopo = 20;

% frequência de amostragem 
fn = fa/2;

% Obter o valor desse angulo predistorcido para compensar a distorção na 
% frequência causada pela transformação bilinear
Ts1 = fs1/fn;
Ts2 = fs2/fn;
Tp1 = fp1/fn;
Tp2 = fp2/fn;

% transformação bilinear 
Ls1 = 2*tan((Ts1*pi)/2);
Ls2 = 2*tan((Ts2*pi)/2);
Lp1 = 2*tan((Tp1*pi)/2);
Lp2 = 2*tan((Tp2*pi)/2);

B = (Lp2 - Lp1);
L0 = sqrt(Lp2*Lp1);

Wp = 1;

% Cálculo do protótipo com Omega
Os1 = abs(B*Ls1 / (Ls1^2 - L0^2));
Os2 = abs(B*Ls2 / (Ls2^2 - L0^2));
Ws = min(Os1,Os2);

% Encontrando a ordem e a freq de corte 
% Determinando Hp(p) prototipo
[n,Wn] = cheb2ord(Wp,Ws,Ap,As,'s');
[bp,ap] = cheby2(n,As, Wn,'s');
figure(1);
subplot(211);
freqs(bp,ap)
title('Protótipo do filtro')

syms p s z;

% Obtendo numerador e denominador
Np(p) = poly2sym(bp, p);
Dp(p) = poly2sym(ap, p);
Hp(p) = Np/Dp;
pretty(vpa(Hp(p),3));

subplot(212);
zplane(bp, ap);
title('Polos e zeros do prototipo');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtro analogico
% Substituição para obter um rejeita-faixa
Hs(s) = subs(Hp,((B*s)/(s^2 + L0^2)));

% Obtendo numerador e denominador
[Ns, Ds] = numden(Hs(s));
bs = sym2poly(Ns);
as = sym2poly(Ds);

% Normalizando
Nbs = bs/as(1);
Nas = as/as(1);

figure(2);
subplot(211);
[h2, w2] = freqs(Nbs, Nas, 10000);
semilogx(w2,20*log10(abs(h2)));
title('Filtro analógico')

subplot(212);
zplane(bs, as);
title('Polos e zeros filtro analógico');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtro Digital
% Substituição para obter um rejeita-faixa
Hz(z)= subs(Hs,(2*((z-1)/(z+1))));
pretty(vpa(Hz(z), 2))

% Obtendo numerador e denominador
[Nz, Dz] = numden(Hz(z));
bz = sym2poly(Nz)*10^(Gtopo/20);
az = sym2poly(Dz);

% Normalizando
Nbz = bz/az(1);
Naz = az/az(1);

figure(3);
subplot(211);
[h3, w3] = freqz(Nbz, Naz);
semilogx(w3/pi*(fa/2),20*log10(abs(h3)));
grid on; 
hold on;

plot(fs1,Gtopo-As,'or')
plot(fp1,Gtopo-Ap,'or')
plot(fp2,Gtopo-Ap,'or')
plot(fs2,Gtopo-As,'or')

line([500 fs1 fs1 fs2 fs2 4000],[Gtopo Gtopo Gtopo-As Gtopo-As Gtopo Gtopo],'Color','red','LineStyle','--');
line([500 fp1 fp1],[Gtopo-Ap Gtopo-Ap -70],'Color','red','LineStyle','--');
line([fp2 fp2 4000],[-70 Gtopo-Ap Gtopo-Ap],'Color','red','LineStyle','--');

axis([500 4000 -70 30])
title('Filtro digital')
xlabel('x');
ylabel('y');

subplot(212);
zplane(bz, az);
title('Polos e zeros filtro digital');
hold off;