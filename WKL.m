clc
clear all

%% Equação da onda - Warming-Kutler-Lomax

%Variáveis
c = 1.0;          %velocidade da onda
L = 40.0;         %comprimento do domínio
NX = 41;          %número de pontos da malha
dx = L / (NX-1);  %passo de espaço 
t = 0.0;          %tempo
t_f = 18.0;       %tempo final
n = 3;            %número da onda
ni = 1;         %número de Courant
dt = ni*(dx/c);   %passo de tempo
k_m = (n*pi)/L;   %número de onda
A_0 = 1;          %amplitude inicial

%Discretização do domínio
i = 1 : NX;
x = (i - 1)*dx;   %domínio espacial
NT = t_f / dt;    %número do loop do tempo

%Condição inicial (para o tempo n)
u_0 = sin(2 * n * pi * ( x / L ));
u_num = u_0;

for k = 1 : NT     %Loop no tempo
    
    %Para o time step 1 (tempo n^(1))
    for j = 1 : NX-1
        u_1(j) = u_0(j)-(2/3)*ni*(u_0(j+1)-u_0(j));
    end
    u_1(NX) = u_1(1);   %condição de contorno periódica
    
    %Para o time step 2 (tempo n^(2))
    for j = 2 : NX
        u_2(j) = (u_0(j)+u_1(j)-(2/3)*ni*(u_1(j)-u_1(j-1)))/2;
    end
    u_2(1) = u_2(NX);
    
    %Para o time step 3 (tempo n+1)
    w = 3;
    for j = 3 : NX-2
        u_num(j) = u_0(j)-(ni*(-2*u_0(j+2)+7*u_0(j+1)-7*u_0(j-1)+2*u_0(j-2)))/24 - (3/8)*ni*(u_2(j+1)-u_2(j-1)) - (w/24)*(u_0(j+2)-4*u_0(j+1)+6*u_0(j)-4*u_0(j-1)+u_0(j-2));
    end
    u_num(1) = u_0(1)-(ni*(-2*u_0(3)+7*u_0(2)-7*u_0(NX-1)+2*u_0(NX-2)))/24 - (3/8)*ni*(u_2(2)-u_2(NX-1)) - (w/24)*(u_0(3)-4*u_0(2)+6*u_0(1)-4*u_0(NX-1)+u_0(NX-2));
    u_num(2) = u_0(2)-(ni*(-2*u_0(4)+7*u_0(3)-7*u_0(1)+2*u_0(NX-1)))/24 - (3/8)*ni*(u_2(3)-u_2(1)) - (w/24)*(u_0(4)-4*u_0(3)+6*u_0(2)-4*u_0(1)+u_0(NX-1));
    u_num(NX-1) = u_0(NX-1)-(ni*(-2*u_0(2)+7*u_0(NX)-7*u_0(NX-2)+2*u_0(NX-3)))/24 - (3/8)*ni*(u_2(NX)-u_2(NX-2)) - (w/24)*(u_0(2)-4*u_0(NX)+6*u_0(NX-1)-4*u_0(NX-2)+u_0(NX-3));
    u_num(NX) = u_num(1);
    
    %atualizando u e t
    t = t + dt;
    u_0 = u_num;
    
    %solução analítica
    u_ana = sin(2 * n * pi * ( (x - c*t) / L ));
    
    %plotando os resultados
    plot(x,u_num,'-bo');
    hold on
    plot(x,u_ana,'-r');
    hold off
    title('Warming-Kutler-Lomax: n=3 e ni=1');
    grid on;
    xlabel('Número de nós');
    ylabel('u(x,t)');
end