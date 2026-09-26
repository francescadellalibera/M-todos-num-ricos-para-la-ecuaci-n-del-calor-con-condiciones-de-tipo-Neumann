clc
clear all
close all

%% =========================
%  Parámetros físicos
%% =========================
save_times = [1 5 10 15 20 25 50 100 200 300]; % en segundos
L = 0.01;              % profundidad piel (m)
Nx = 100;
dx = L/(Nx-1);

kth = 0.37;            % W/(m·K)
rho = 1000;            % kg/m^3
c = 3500;              % J/(kg·K)

alpha = kth/(rho*c);   % difusividad térmica (m^2/s)

%% =========================
%  Flujo en superficie
%% =========================
g0 = 2000;      % W/m^2
kappa = 0.05;   % 1/s

%% =========================
%  Condición inicial
%% =========================
T_init = 310;   % K (37°C)
T = T_init * ones(Nx,1);
T_new = T;

%% =========================
%  Tiempo
%% =========================
%dt = 0.2 * dx^2 / alpha;% estabilidad Euler explícito
dt = 0.05 * dx^2 / alpha;
T_end = 300;               % s
Nt = round(T_end/dt);

x = linspace(0,L,Nx);

%% =========================
%  Bucle temporal
%% =========================
for n = 1:Nt
    
    t = (n-1)*dt;
    
    % flujo superficial
    %g = g0 * exp(-kappa * t);
    g = g0 * (1 - exp(-kappa*t)) * exp(-kappa*t);
    
    %% -------------------------
    % Interior (difusión)
    %% -------------------------
    for i = 2:Nx-1
        T_new(i) = T(i) + alpha*dt/dx^2 * ...
            (T(i+1) - 2*T(i) + T(i-1));
    end
    
    %% -------------------------
    % Frontera x = 0 (Neumann con flujo)
    %% -------------------------
    T_new(1) = T(1) ...
        + 2*alpha*dt/dx^2 * (T(2) - T(1)) ...
        + 2*dt/(rho*c*dx) * g;
    
    %% -------------------------
    % Frontera x = L (aislado)
    %% -------------------------
    T_new(Nx) = T_new(Nx-1);
    
    %% actualizar
    T = T_new;
    
    %% -------------------------
    % visualización
    %% -------------------------
    if mod(n,50)==0
        plot(x, T-273, 'LineWidth', 2);
        xlabel('Profundidad (m)');
        ylabel('Temperatura (°C)');
        title(['t = ', num2str(t,'%.1f'), ' s']);
        %ylim([20 40]);
       ylim([20 max(50, max(T-273)+1)]);
        grid on;
        drawnow;
    end
     % Guardar frames en tiempos específicos
   for k = 1:length(save_times)
        if abs(t - save_times(k)) < dt/2
            figure(1)
            plot(x, T-273, 'LineWidth', 2);
            xlabel('Profundidad (m)');
          ylabel('Temperatura (°C)');
            title(['t = ', num2str(t,'%.1f'), ' s']);
            ylim([20 50]);
            xlim([0 L]);
            grid on;
            
            filename = ['frame_t_', num2str(round(t)), '.png'];
            saveas(gcf, filename);
        end
   end
end

