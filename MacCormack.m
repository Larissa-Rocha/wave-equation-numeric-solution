clc
clear all
%% Equa��o da onda - MacCormack

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

%Condi��o inicial (Para t=0, n�vel n)
u_0 = sin(2 * n * pi * ( x / L ));
u = u_0;
u_num = u_0;

%Solu��o
%Loop no tempo
NT = t_f / dt;
for k = 1 : NT
    
    %c�lculo do preditor
    for j = 1 : NX
        if j <= NX-1
            u_pred(j) = u_0(j) - (ni)*(u_0(j+1)-u_0(j));
        else
            u_pred(j) = u_0(j) - (ni)*(u_0(2)-u_0(j));       %condi��o de contorno peri�dica          
        end
    end
    u_0 = u_pred;
    
    %c�lculo do corretor
    for j = 1 : NX
        if j == 1
            u_num(j) = (u(j) + u_pred(j) - (ni)*(u_pred(j)-u_pred(NX-1)))/2;    %condi��o de contorno peri�dica
        else
            u_num(j) = (u(j) + u_pred(j) - (ni)*(u_pred(j)-u_pred(j-1)))/2;
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
    plot(u, '-bo');
    hold off
    title('MacCormack: n=3 e ni=1');
    grid on;
    xlabel('N�mero de n�s');
    ylabel('u(x,t)');
end