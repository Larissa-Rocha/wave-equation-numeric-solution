clc
clear all

%% Equa��o da onda - Leap frog

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
i = 1 : NX;
x = (i - 1)*dx;   %dom�nio espacial
NT = t_f / dt;    %n�mero de tempo final

u_0 = zeros(1,NX);
u_store = zeros(NT,NX);
%Condi��o inicial (para o tempo n-1)
u_0 = sin(2 * n * pi * ( x / L ));
u = u_0;
u_store(1,:) = u_0;   %armazenando o primeiro tempo

% para o primeiro time step (n)
% usando FTCS e condi��es de contorno peri�dicas
for j = 1 : NX
    if j == 1
        u(1,j) = u_0(1,j)-(ni/2)*(u_0(1,j+1) - u_0(1,NX-1));  %contorno peri�dico
    elseif j <= NX-1
        u(1,j) = u_0(1,j)-(ni/2)*(u_0(1,j+1)-u_0(1,j-1));
    else 
        u(1,j) = u_0(1,j)-(ni/2)*(u_0(1,2)-u_0(1,j-1));     %contorno peri�dico
    end

end
u_0 = u;
u_store(2,:) = u;

%para o segundo time step (n+1) em diante
%usando Leap Frog e condi��es de contorno peri�dicas

for k = 3 : NT     %Loop no tempo
    for j = 1 : NX
        if j == 1
            u(1,j)=u_store(k-2,j)-ni*(u_0(1,j+1)-u_0(1,NX-1));    %contorno peri�dico
        elseif j<= NX-1
            u(1,j)=u_store(k-2,j)-ni*(u_0(1,j+1)-u_0(1,j-1));
        else
            u(1,j)=u_store(k-2,j)-ni*(u_0(1,2)-u_0(1,j-1));     %contorno peri�dico
        end
    end
    t = t + dt;
    u_0 = u;
    u_store(k,:) = u;
    
    %solu��o anal�tica
    u_ana = sin(2 * n * pi * ( (x - c*t) / L ));
    
    %plotando os resultados
    plot(x,u,'-bo');
    hold on
    plot(x,u_ana,'-r');
    hold off
    title('Leap Frog: n=3 e ni=1');
    grid on;
    xlabel('N�mero de n�s');
    ylabel('u(x,t)');
end

%%
%Determinando beta
beta = k_m * dx;

%Determinando os erros na amplitude e na fase
G = 1;                                              %M�dulo do Fator de amplifica��o
G_n = (1 - G^t_f)*A_0;                              %ap�s 18 passos de tempo
Fi = atan2(-ni*sin(beta), 1 - ni + ni*cos(beta));   %�ngulo de fase
Fi_e = -beta*ni;                                    %�ngulo de fase exato
Fi_n = t_f*(Fi_e - Fi);                             %ap�s 18 passos de tempo