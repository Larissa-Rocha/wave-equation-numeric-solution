clc
clear all
%% Wave equation - Lax-Wendroff

%Variables
c = 1.0;          %wave velocity
L = 40.0;         %lenght
NX = 41;          %number of nodes
dx = L / (NX-1);  %space step 
t = 0.0;          %time
t_f = 18.0;       %final time
n = 3;            %wave number
ni = 1;         %Courant number
dt = ni*(dx/c);   %time step
k_m = (n*pi)/L;   %wave number
A_0 = 1;          %initial amplitude

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
    
    for j = 2 : NX-1    %Loop no espaço
        u_num(j) = u(j)-(ni/2)*(u(j+1)-u(j-1))-(ni^2/2)*(u(j+1)-2*u(j)+u(j-1)); %Esquema Lax-Wendroff  
    end
    
    %Condições de contorno periódicas
    u_num(1) = u(1)-(ni/2)*(u(2)-u(NX-1))-(ni^2/2)*(u(2)-2*u(1)+u(NX-1));
    u_num(NX) = u_num(1);
    
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
    title('Lax Wendroff: n=3 e ni=1');
    grid on;
    xlabel('Número de nós');
    ylabel('u(x,t)');
end
