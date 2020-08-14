clc
clear all
%% Equação da onda - Upwind de segunda ordem

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

%Condição inicial (Para t=0)
u_0 = sin(2 * n * pi * ( x / L ));
u = u_0;
u_num = u_0;

%Solução
%Loop no tempo
NT = t_f / dt;
for k = 1 : NT
    for j = 2 : NX
        if j==2
            u_num(j) = u(j) - ((ni) * (u(j) - u(j-1))) + (ni * (ni-1) * (u(j)-2*u(j-1)+u(NX-1)))/2;
        else
            u_num(j) = u(j) - ((ni) * (u(j) - u(j-1))) + (ni * (ni-1) * (u(j)-2*u(j-1)+u(j-2)))/2;
        end
    end
    u_num(1) = u_num(NX);   %condição de contorno periódica
    
    %Atualizar t e u
    t = t + dt;
    u = u_num;
    
    %Solução analítica
    u_ana = sin ( 2 * n * pi * ((x - c*t) / L));
    
    %Plotando os resultados
      plot(u_ana, 'r-');
      hold on
      plot(u_num, '-bo');
      hold off
      title('Upwind de 2ª ordem: n=3 e ni=1');
    grid on;
    xlabel('Número de nós');
    ylabel('u(x,t)');
end