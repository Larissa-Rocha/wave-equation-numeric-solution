clc
clear all
%% Equação da onda - Upstream

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
x_min = 1;
x_max = NX-1;
x = x_min - dx : dx : x_max;     %domínio espacial

%Condição inicial (Para t=0)
u_0 = sin(2 * n * pi * ( x / L ));
u = u_0;
u_num = u_0;

%Solução
%Loop no tempo
NT = t_f / dt;
for k = 1 : NT
    
    %Condições de contorno periódicas
    %u(NX-1) = u(0);
%     u(1) = u(NX);
    
    %Esquema Upstream
    %Loop no espaço
    for j = 1 : NX
        if (j == 1)
            u_num(j) = u(j)*(1 - ni) + ni*u(NX-1);
        else
            u_num(j) = u(j)*(1 - ni) + ni*u(j-1);
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
    plot(u, 'bo');
    hold off
    title('Upstream: n=3 e ni=1');
    grid on;
    xlabel('Número de nós');
    ylabel('u(x,t)');
end

%%
%Determinando beta
beta = k_m * dx;

%Determinando os erros na amplitude e na fase
G = sqrt((1 - ni + ni*cos(beta))^2 + (-ni*sin(beta))^2);    %Fator de amplificação
G_n = (1 - G^t_f)*A_0;                                      %após 18 passos de tempo
Fi = atan2(-ni*sin(beta), 1 - ni + ni*cos(beta));           %ângulo de fase
Fi_e = -beta*ni;                                            %ângulo de fase exato
Fi_n = t_f*(Fi_e - Fi);                                     %após 18 passos de tempo

