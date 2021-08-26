% Filtro BS (rejeita-faixa) 
close all;
clear all;
clc;

% Dados do filtro
fp1 = 1200; % é a frequência de passagem em Hz
fs1 = 2000; % é a frequência de rejeição em Hz (stopband)
fs2 = 3000; % é a frequência de rejeição em Hz (stopband)
fp2 = 3800; % é a frequência de passagem em Hz
fa = 8000;  % é a frequência de amostragem
Ap = 3;     % é a atenuação em dB na frequência de passagem      
As = 40;    % é a atenuação em dB na frequência de stopband   
Gtopo = 20;

Os1 = fs1/(fa/2);
Op1 = fp1/(fa/2);
Op2 = fp2/(fa/2);
Os2 = fs2/(fa/2);

% Parks McClellan
f = [fp1 fs1 fs2 fp2];
a = [1 0 1];
desvio = [(10^((Ap/2)/20)-1)/(10^((Ap/2)/20)+1) 10^(-As/20) (10^((Ap/2)/20)-1)/(10^((Ap/2)/20)+1)];
[N1,fo,ao,w1] = firpmord(f,a,desvio,fa);

% Ajustes no ganho
b = firpm(N1,fo,ao,w1);
Glin1 = 10^((Gtopo-1.15)/20);
bAjustado = b*Glin1;

figure(1) 
[h1,w1] = freqz(bAjustado,1,1024,fa);
plot(w1, 20*log10(abs(h1)));
grid on;
title(sprintf('Resposta em Frequência n = %d - Parks McClellan',N1))
ylabel('dB')
xlabel('Frequencia (Hz)')
axis([0 fa/2 -As-15 Gtopo+5]);

hold on;

line([0 fs1],[Gtopo Gtopo],'Color','magenta','LineStyle','--')
line([0 fp1],[Gtopo-Ap Gtopo-Ap],'Color','magenta','LineStyle','--')
line([fp1 fp1],[-As-15 Gtopo-Ap],'Color','magenta','LineStyle','--')
line([fs1 fs1],[-As+Gtopo Gtopo],'Color','magenta','LineStyle','--')

line([fs1 fs2],[-As+Gtopo -As+Gtopo],'Color','magenta','LineStyle','--')
line([fp2 fp2],[-As-15 Gtopo-Ap],'Color','magenta','LineStyle','--')
line([fs2 fs2],[-As+Gtopo Gtopo],'Color','magenta','LineStyle','--')
line([fs2 fa/2],[Gtopo Gtopo],'Color','magenta','LineStyle','--')
line([fp2 fa/2],[Gtopo-Ap Gtopo-Ap],'Color','magenta','LineStyle','--')

figure(2)
grpdelay(bAjustado);
title('Atraso de grupo - Parks McClellan');

figure(3);
zplane(bAjustado);
title('Pólos e Zeros - Parks McClellan');
xlabel('Re');
ylabel('Im');
grid on;

figure(4);
stem(bAjustado);
title('Resposta ao impulso - Parks McClellan');
grid on;

%Janela Fixa(hamming)
N2 = 30;
M = N2/2;

ws1 = Os1*pi; 
ws2 = Os2*pi;
wp1 = Op1*pi; 
wp2 = Op2*pi;

wc1 = (ws1+wp1)/2;
wc2 = (ws2+wp2)/2;

n = [-M:M];
CBS = (sin(wc1*n) - sin(wc2*n))./(pi*n);
CBS(M+1) = 1 - (wc2-wc1)/pi;

FiltroHamming = CBS.*hamming(2*M+1)';
FiltroHammingAjustado = FiltroHamming*10^((Gtopo)/20);

figure(5);
[h2, w2] = freqz(FiltroHammingAjustado,1);
plot((w2/pi)*(fa/2),20*log10(abs(h2)));
grid on;
hold on;

line([0 fs1],[Gtopo Gtopo],'Color','magenta','LineStyle','--')
line([0 fp1],[Gtopo-Ap Gtopo-Ap],'Color','magenta','LineStyle','--')
line([fp1 fp1],[-As-60 Gtopo-Ap],'Color','magenta','LineStyle','--')
line([fs1 fs1],[-As+Gtopo Gtopo],'Color','magenta','LineStyle','--')

line([fs1 fs2],[-As+Gtopo -As+Gtopo],'Color','magenta','LineStyle','--')
line([fp2 fp2],[-As-60 Gtopo-Ap],'Color','magenta','LineStyle','--')
line([fs2 fs2],[-As+Gtopo Gtopo],'Color','magenta','LineStyle','--')
line([fs2 fa/2],[Gtopo Gtopo],'Color','magenta','LineStyle','--')
line([fp2 fa/2],[Gtopo-Ap Gtopo-Ap],'Color','magenta','LineStyle','--')
hold off;

title(sprintf('Resposta em Frequência n = %d - Hamming',N2))
axis([0 fa/2 -As-60 Gtopo+5]);


figure(6)
grpdelay(FiltroHammingAjustado);
title('Atraso de grupo - Hamming');

figure(7);
zplane(FiltroHammingAjustado);
title('Pólos e Zeros - Hamming');
xlabel('Re');
ylabel('Im');
grid on;

figure(8);
stem(FiltroHammingAjustado);
title('Resposta ao impulso - Hamming');
grid on;

%%
% Janela Fixa(hamming) com equação
N3 = 30;

wc1 = (Os1+Op1)/2;
wc2 = (Os2+Op2)/2;


Wn = [wc1 wc2];

h_fir = fir1(N3,Wn,'stop',hamming(N3+1));

% Ajustes do ganho
Glin2 = 10^((Gtopo-2)/20);
h_firAjustado = h_fir*Glin2;
[hh,wh] = freqz(h_firAjustado);
figure(9);
plot(wh*fa/2/pi,20*log10(abs(hh)));
axis([0 fa/2 -As-20 Gtopo+5]);
hold on;

line([0 fs1],[Gtopo Gtopo],'Color','magenta','LineStyle','--')
line([0 fp1],[Gtopo-Ap Gtopo-Ap],'Color','magenta','LineStyle','--')
line([fp1 fp1],[-As-25 Gtopo-Ap],'Color','magenta','LineStyle','--')
line([fs1 fs1],[-As+Gtopo Gtopo],'Color','magenta','LineStyle','--')

line([fs1 fs2],[-As+Gtopo -As+Gtopo],'Color','magenta','LineStyle','--')
line([fp2 fp2],[-As-25 Gtopo-Ap],'Color','magenta','LineStyle','--')
line([fs2 fs2],[-As+Gtopo Gtopo],'Color','magenta','LineStyle','--')
line([fs2 fa/2],[Gtopo Gtopo],'Color','magenta','LineStyle','--')
line([fp2 fa/2],[Gtopo-Ap Gtopo-Ap],'Color','magenta','LineStyle','--')
hold off;
title(['Janela fixa hamming n = ' num2str(N2)])

%% TODOS OS FILTROS JUNTOS

figure(10);
plot(w1, 20*log10(abs(h1)));
hold on;
plot(wh*fa/2/pi,20*log10(abs(hh)));
hold on;
plot((w2/pi)*(fa/2),20*log10(abs(h2)));
legend(['Parks McClellan n = ', num2str(N1)],['Janela fixa hamming n = ', num2str(N2)],['Janela fixa hamming - com equação n = ', num2str(N3)]);

axis([0 fa/2 -As-20 Gtopo+5]);
line([0 fs1],[Gtopo Gtopo],'Color','magenta','LineStyle','--')
line([0 fp1],[Gtopo-Ap Gtopo-Ap],'Color','magenta','LineStyle','--')
line([fp1 fp1],[-As-25 Gtopo-Ap],'Color','magenta','LineStyle','--')
line([fs1 fs1],[-As+Gtopo Gtopo],'Color','magenta','LineStyle','--')

line([fs1 fs2],[-As+Gtopo -As+Gtopo],'Color','magenta','LineStyle','--')
line([fp2 fp2],[-As-25 Gtopo-Ap],'Color','magenta','LineStyle','--')
line([fs2 fs2],[-As+Gtopo Gtopo],'Color','magenta','LineStyle','--')
line([fs2 fa/2],[Gtopo Gtopo],'Color','magenta','LineStyle','--')
line([fp2 fa/2],[Gtopo-Ap Gtopo-Ap],'Color','magenta','LineStyle','--')
hold off;



