clc
clear all
%% Equa��o da onda - Upstream

%Vari�veis
c = 1.0;          %velocidade da onda
L = 40.0;         %comprimento do dom�nio
NX = 41;          %n�mero de pontos da malha
dx = L / (NX-1);  %passo de espa�o 
t = 0.0;          %tempo
t_f = 18.0;       %tempo final
n = 3;            %n�mero da onda
ni = 1;         %n�mero de Courant
dt = ni*(dx/c);   %passo de tempo
k_m = (n*pi)/L;   %n�mero de onda
A_0 = 1;          %amplitude inicial

%Discretiza��o do dom�nio
x_min = 1;
x_max = NX-1;
x = x_min - dx : dx : x_max;     %dom�nio espacial

%Condi��o inicial (Para t=0)
u_0 = sin(2 * n * pi * ( x / L ));
u = u_0;
u_num = u_0;

%Solu��o
%Loop no tempo
NT = t_f / dt;
for k = 1 : NT
    
    %Condi��es de contorno peri�dicas
    %u(NX-1) = u(0);
%     u(1) = u(NX);
    
    %Esquema Upstream
    %Loop no espa�o
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
    
    %Solu��o anal�tica
    u_ana = sin ( 2 * n * pi * ((x - c*t) / L));
    
    %Plotando os resultados
    plot(u_ana, 'r-');
    hold on
    plot(u, 'bo');
    hold off
    title('Upstream: n=3 e ni=1');
    grid on;
    xlabel('N�mero de n�s');
    ylabel('u(x,t)');
end

%%
%Determinando beta
beta = k_m * dx;

%Determinando os erros na amplitude e na fase
G = sqrt((1 - ni + ni*cos(beta))^2 + (-ni*sin(beta))^2);    %Fator de amplifica��o
G_n = (1 - G^t_f)*A_0;                                      %ap�s 18 passos de tempo
Fi = atan2(-ni*sin(beta), 1 - ni + ni*cos(beta));           %�ngulo de fase
Fi_e = -beta*ni;                                            %�ngulo de fase exato
Fi_n = t_f*(Fi_e - Fi);                                     %ap�s 18 passos de tempo

