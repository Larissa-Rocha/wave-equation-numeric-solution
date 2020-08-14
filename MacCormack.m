clc
clear all
%% Equação da onda - MacCormack

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

%Condição inicial (Para t=0, nível n)
u_0 = sin(2 * n * pi * ( x / L ));
u = u_0;
u_num = u_0;

%Solução
%Loop no tempo
NT = t_f / dt;
for k = 1 : NT
    
    %cálculo do preditor
    for j = 1 : NX
        if j <= NX-1
            u_pred(j) = u_0(j) - (ni)*(u_0(j+1)-u_0(j));
        else
            u_pred(j) = u_0(j) - (ni)*(u_0(2)-u_0(j));       %condição de contorno periódica          
        end
    end
    u_0 = u_pred;
    
    %cálculo do corretor
    for j = 1 : NX
        if j == 1
            u_num(j) = (u(j) + u_pred(j) - (ni)*(u_pred(j)-u_pred(NX-1)))/2;    %condição de contorno periódica
        else
            u_num(j) = (u(j) + u_pred(j) - (ni)*(u_pred(j)-u_pred(j-1)))/2;
        end
    end
    
    %Atualizar t e u
    t = t + dt;
    u = u_num;
    
    %Solução analítica
    u_ana = sin ( 2 * n * pi * ((x - c*t) / L));
    
    %Plotando os resultados
    plot(u_ana, 'r-');
    hold on
    plot(u, '-bo');
    hold off
    title('MacCormack: n=3 e ni=1');
    grid on;
    xlabel('Número de nós');
    ylabel('u(x,t)');
end