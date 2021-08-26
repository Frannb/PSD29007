% Filtro HP (Passa-alta)
close all;
clear all;
clc;

% Dados do filtro
fs = 40;        % é a frequência de rejeição em Hz (stopband)
fp = 200;       % é a frequência de passagem em Hz
Ap = 1;         % é a atenuação em dB na frequência de passagem
As = 50;        % é a atenuação em dB na frequência de stopband
fa = 3*fp;      % é a frequência de amostragem
Gtopo = 10;

%Parks-McClellan
f = [fs fp];    
a = [0 1];        
desvio = [10^(-As/20) (10^(Ap/20)-1)/(10^(Ap/20)+1)]; 
[n1,fo,ao,w1] = firpmord(f,a,desvio,fa);
b = firpm(n1,fo,ao,w1);

% Ajustes no ganho
Glin1 = 10^((Gtopo-0.5)/20);
bAjustado = b*Glin1;

[h1,w1] = freqz(bAjustado,1,1024,fa);

figure(1);
plot((w1), 20*log10(abs(h1)));
title(['Filtro Parks McClellan n = ' num2str(n1)])
ylabel('dB')
xlabel('Frequencia (Hz)')
axis([0 fa/2 -As-10 Gtopo+5]);
hold on;

line([0 fs],[-As+Gtopo -As+Gtopo],'Color','magenta','LineStyle','--')
line([fs fs],[-As+Gtopo Gtopo],'Color','magenta','LineStyle','--')
line([fs fa/2],[Gtopo Gtopo],'Color','magenta','LineStyle','--')
line([fp fa/2],[Gtopo-Ap Gtopo-Ap],'Color','magenta','LineStyle','--')
line([fp fp],[-As-10 Gtopo-Ap],'Color','magenta','LineStyle','--')

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


% Janela Ajustavel
[n2,wn,beta,ftype] = kaiserord(f,a,desvio,fa);

wkaiser = kaiser(n2+1,beta);
h_fir = fir1(n2,wn,ftype,wkaiser,'noscale');

% Ajustes no ganho
Glin2 = 10^((Gtopo - 0.5)/20);
h_firAjustado = h_fir*Glin2;

[h2,w2] = freqz(h_firAjustado);

figure(5);
plot(w2*fa/2/pi,20*log10(abs(h2)));
ylabel('dB')
xlabel('Frequencia (Hz)')
title(['Filtro kaiser n = ' num2str(n2)]);
axis([0 300 -85 11]);
hold on;

line([0 fs],[-As+Gtopo -As+Gtopo],'Color','magenta','LineStyle','--')
line([fs fs],[-As+Gtopo Gtopo],'Color','magenta','LineStyle','--')
line([fs fa/2],[Gtopo Gtopo],'Color','magenta','LineStyle','--')
line([fp fa/2],[Gtopo-Ap Gtopo-Ap],'Color','magenta','LineStyle','--')
line([fp fp],[-As-30 Gtopo-Ap],'Color','magenta','LineStyle','--')

figure(6)
grpdelay(h_firAjustado);
title('Atraso de grupo - Kaiser');

figure(7);
zplane(h_firAjustado);
title('Pólos e Zeros - Kaiser');
xlabel('Re');
ylabel('Im');
grid on;

figure(8);
stem(h_firAjustado);
title('Resposta ao impulso - Kaiser');
grid on;

figure(9);
plot(w1,20*log10(abs(h1)));
hold on
plot(w2*fa/2/pi,20*log10(abs(h2)));
title(['Filtro Parks McClellan n = ' num2str(n1) ' e Filtro kaiser n = ' num2str(n2)]);
axis([0 300 -85 11]);
line([0 fs],[-As+Gtopo -As+Gtopo],'Color','magenta','LineStyle','--')
line([fs fs],[-As+Gtopo Gtopo],'Color','magenta','LineStyle','--')
line([fs fa/2],[Gtopo Gtopo],'Color','magenta','LineStyle','--')
line([fp fa/2],[Gtopo-Ap Gtopo-Ap],'Color','magenta','LineStyle','--')
line([fp fp],[-As-40 Gtopo-Ap],'Color','magenta','LineStyle','--')
ylabel('dB')
xlabel('Frequencia (Hz)')
hold off;