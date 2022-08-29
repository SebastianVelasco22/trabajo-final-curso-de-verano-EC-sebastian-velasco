clear all;
close all;
clc;

%% inicializar variables:
fs = 150;
t = linspace(-5,5,fs);
ax = [-4 4 -1 5];
ax_tau = [-10 10 -1 5];
ax_conv = [-10 10 -1 12];

%% seniales f(t) y g(t):
disp('Función o señal: f(t)')
disp('Función o señal: g(t')
f = ((3/2)*t+3).*(t>= -2 & t<=0) + (-(3/2)*t+3).*(t> 0 & t<=2); %(t>=0 & t<=3).*(t.^2);
g = 4*exp(-t).*(t>=0); %(t+1).*(t>= -1 & t<=0) + (-t+1).*(t> 0 & t<=1); 

figure(1)
plot(t,f,'Color','blue','LineWidth',2),axis(ax);
grid on;
xlabel('\tau','FontSize', 16)
ylabel('f(\tau)','FontSize', 16)

% g=fliplr(g); %invierte  la señal: g(-t)% no se uso realmente
figure(2)
plot(t,g,'Color','red','LineWidth',2),axis(ax);
grid on;
xlabel('\tau','FontSize', 16)
ylabel('g(\tau)','FontSize', 16)

disp('Presiona una tecla ...');
pause;

%% Cambio de t a Tau para convolucion:
tau = linspace(min(t)*2, max(t)*2,length(f) + length(g)+ 1);
f_tau = ((3/2)*tau+3).*(tau>= -2 & tau<=0) + (-(3/2)*tau+3).*(tau> 0 & tau<=2);
g_tau = 4*exp(-tau).*(tau>=0);
iter = length(zeros(1,length(f) + length(g)+ 1));
convolution = zeros(1,length(f) + length(g)+ 1);

%% Convolucion a pasos de f con g:
figure(3)
for i = 1:iter 
    mov = min(tau) + (i*(tau(2) - tau(1))); %(10*i-fs)/fs;
    g_modif = 4*exp(-(mov-tau)).*(mov-tau>=0);
    convolution(i) = trapz(tau, f_tau.*g_modif);%Hace integral de funciones
   
    subplot(2,2,1:2)
    hold off;
    plot(tau,f_tau,'Color','blue','LineWidth',2),axis(ax_tau);
    hold on; 
    plot(tau,g_modif,'Color','red' ,'LineWidth',2),axis(ax_tau);
    grid on;
    xlabel('\tau', 'FontSize', 16);
    ylabel('g(t - \tau)', 'FontSize', 16);
   
    subplot(2,2,3:4)
    hold off
    plot(tau(1:i),convolution(1:i),'Color','black','LineWidth',2),axis(ax_conv);
    grid on;
    xlabel('\tau', 'FontSize', 16); 
    ylabel('\int f(\tau)g(t - \tau)d\tau', 'FontSize', 16);
    
    pause(0.01)
end
%% Convolucion a pasos de g con f:
figure(4)
for i = 1:iter 
    mov = min(tau) + (i*(tau(2) - tau(1))); %(10*i-fs)/fs;
    f_modif = ((3/2)*(mov-tau)+3).*((mov-tau)>= -2 & (mov-tau)<=0) + (-(3/2)*(mov-tau)+3).*((mov-tau)> 0 & (mov-tau)<=2);
    convolution(i) = trapz(tau, g_tau.*f_modif);%Hace integral de funciones
   
    subplot(2,2,1:2)
    hold off;
    plot(tau,g_tau,'Color','red','LineWidth',2),axis(ax_tau);
    hold on; 
    plot(tau,f_modif,'Color','blue' ,'LineWidth',2),axis(ax_tau);
    grid on;
    xlabel('\tau', 'FontSize', 16);
    ylabel('f(t - \tau)', 'FontSize', 16);
   
    subplot(2,2,3:4)
    hold off
    plot(tau(1:i),convolution(1:i),'Color','black','LineWidth',2),axis(ax_conv);
    grid on;
    xlabel('\tau', 'FontSize', 16); 
    ylabel('\int g(\tau)f(t - \tau)d\tau', 'FontSize', 16);
    
    pause(0.01)
end
%% Convolucion:
figure(5)
plot(tau,convolution,'Color','black','LineWidth',2),axis(ax_conv);
grid on;
xlabel('t', 'FontSize', 16); 
ylabel('f(t) * g(t)', 'FontSize', 16);

%% Titulos de frames:
figure(1)
title('Señal o Función f(\tau) ');
figure(2)
title('Señal o Función g(\tau) ');
figure(5)
title('Convolución ');

