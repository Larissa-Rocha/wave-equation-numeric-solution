clc
clear all
%% Equa��o da onda - Lax

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

%Condi��o inicial (Para t=0)
u_0 = sin(2 * n * pi * ( x / L ));
u = u_0;
u_num = u_0;

%Solu��o
%Loop no tempo
NT = t_f / dt;
for k = 1 : NT
    
    %condi��es de contorno peri�dicas
    u_num(1) = ((1-ni)*u(2) + (1+ni)*u(NX-1)) / 2;
    u_num(NX) = u_num(1);
    
    %Esquema Lax
    %Loop no espa�o
    for  j = 2 : NX - 1
        u_num(j) = ((1-ni)*u(j+1) + (1+ni)*u(j-1)) / 2;
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
      title('Lax: n=3 e ni=1');
    grid on;
    xlabel('N�mero de n�s');
    ylabel('u(x,t)');
end

%% 
%Determinando beta
beta = k_m * dx;

%Determinando os erros na amplitude e na fase
G = sqrt((cos(beta))^2 + ni^2 * (sin(beta))^2);             %Fator de amplifica��o
G_n = (1 - G^t_f)*A_0;                                      %ap�s 18 passos de tempo
Fi = atan2(-ni*sin(beta), cos(beta));                       %�ngulo de fase
Fi_e = -beta*ni;                                            %�ngulo de fase exato
Fi_n = t_f*(Fi_e - Fi);                                     %ap�s 18 passos de tempo
