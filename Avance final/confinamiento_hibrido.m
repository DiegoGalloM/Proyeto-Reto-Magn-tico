% confinamiento_hibrido.m
%
% Confinamiento electromagnetico hibrido: F = q(E + v x B)
%
% Etapa 1 — Coulombiana  : N cargas fijas crean barrera de potencial electrico.
% Etapa 2 — Larmor       : Campo magnetico uniforme Bz_0 obliga a giros de Larmor.
% Etapa 3 — Ampere       : Control activo sectorizado — si |r| > R_limit se
%                          activan corrientes de rescate I_high en el cable mas
%                          cercano y sus dos vecinos (parpadeo gris/naranja/rojo).
%
% Hereda:
%   - Geometria y parametros del avance electrico (principal_trampa.m)
%   - Estructura de control magnetico de controlDinamicoMag.m
%
% Integracion: RK4 sobre el estado s = [x, y, vx, vy].
% Conservacion: E_total = Ek + Ep_electrica debe ser ~constante
%               (el campo B no hace trabajo; el campo E si, pero conserva E_mec).

clear; clc; close all

%% 1. PARAMETROS

% --- Particula (electron — consistente con avance electrico) ---
q   =  1.6e-19;      % carga positiva [C] — misma firma que las cargas barrera
m   =  9.11e-31;     % masa de electron [kg]
K   =  8.98755e9;    % constante de Coulomb [N m^2/C^2]
mu0 =  4*pi*1e-7;    % permeabilidad del vacio [T m/A]

% --- Condicion inicial ---
v0_mag = 9e6;                             % rapidez inicial [m/s]
use_fixed_seed = true;
%Para definir una trayectoria en específico
if use_fixed_seed
    rng(1);
else
    rng('shuffle');
end
ang_v = rand * 2*pi;
v0 = [v0_mag*cos(ang_v),  v0_mag*sin(ang_v)];
% r0 alejado del centro para que la particula tenga mas espacio de maniobra
% y los giros de Larmor y el control activo sean visibles en la trayectoria.
% Para comparar con el avance electrico usa r0 = [0.025, 0.02]
r0 = [0.5, 0.0];                          % posicion inicial [m]

% --- Etapa 1: cargas electricas (trampa Coulombiana) ---
N_cargas = 8;
R_cargas = 1.0;                           % radio del anillo de cargas [m]
Q_total  = 8e-6;                          % carga total positiva [C]
Q = ones(1, N_cargas) * Q_total / N_cargas;   % carga individual [C]
ang_Q = 2*pi*(0:N_cargas-1) / N_cargas;
xQ = R_cargas * cos(ang_Q);
yQ = R_cargas * sin(ang_Q);

% --- Etapa 2: campo magnetico uniforme ---
% r_Larmor = m*v / (q*Bz_0)
% Con m=9.11e-31, v=9e6, q=1.6e-19:
%   Bz_0 = 3e-4 T  -> r_L ~ 0.17 m  (giros moderados)
%   Bz_0 = 1e-3 T  -> r_L ~ 0.05 m  (giros apretados, efecto mas visible)
Bz_0 = 1e-3;   % [T]

% --- Etapa 3: cables magneticos (control activo sectorizado) ---
% Los cables comparten el anillo con las cargas electricas.
N_cables = N_cargas;
ang_C = ang_Q;
xC    = xQ;
yC    = yQ;

I_idle   = 10;      % corriente en reposo [A]
I_rescue = 3000;    % corriente de rescate (cable mas cercano) [A]
%                     vecinos reciben 0.5 * I_rescue

R_limit  = R_cargas * 0.50;  % umbral de activacion del rescate [m] = R/2
R_escape = R_cargas;         % radio de escape (toca el anillo) [m]

% --- Tiempo ---
dt    = 1e-10;     % paso de tiempo [s]
steps = 100000;    % pasos totales  -> 10 us de simulacion

% --- Exportar GIF de la animacion ---
% Cambia a true solo cuando quieras guardar el GIF final
make_gif = true;
gif_file = 'confinamiento_hibrido.gif';

% --- Empaquetar todo en estructura para pasar a funciones locales ---
p.q        = q;        p.m        = m;
p.K        = K;        p.mu0      = mu0;
p.N_cargas = N_cargas; p.xQ       = xQ;  p.yQ = yQ; p.Q = Q;
p.Bz_0     = Bz_0;
p.N_cables = N_cables; p.xC       = xC;  p.yC = yC; p.ang_C = ang_C;
p.I_idle   = I_idle;   p.I_rescue = I_rescue;
p.R_limit  = R_limit;  p.R_escape = R_escape;

%% 2. SIMULACION

%Radio larmor teórico
rL_teo = m * v0_mag / (q * Bz_0);
fprintf('\n=== Confinamiento Hibrido Electromagnetico ===\n')
fprintf('Particula : q = %.2e C | m = %.2e kg\n', q, m)
fprintf('v0 = %.2e m/s | dir = %.2f deg\n', v0_mag, rad2deg(ang_v))
fprintf('Bz_0 = %.2e T | r_Larmor teorico = %.4f m\n', Bz_0, rL_teo)
fprintf('R_cargas = %.2f m | R_limit = %.2f m\n\n', R_cargas, R_limit)

[traj, Ek, Ep, t_cpu] = simular(r0, v0, dt, steps, p);

t     = (0:steps-1).' * dt;
r_cil = sqrt(traj(:,1).^2 + traj(:,2).^2);
Etot  = Ek + Ep;
dE    = 100 * (Etot - Etot(1)) / abs(Etot(1));   % deriva relativa [%]

%Condición 1 de escape: El radio del círculo de trayectoria es mayor al de
%escape
escapado  = any(r_cil >= R_escape);
ctrl_frac = 100 * mean(r_cil > R_limit);

fprintf('--- Resultados ---\n')
fprintf('Tiempo CPU       : %.3f s\n', t_cpu)
fprintf('Escapo           : %s\n',    mat2str(escapado))
fprintf('Rescate activo   : %.2f %% del tiempo simulado\n', ctrl_frac)
fprintf('Deriva E_total   : %.5f %% (ideal = 0 %%, B no hace trabajo)\n', dE(end))
fprintf('r_max alcanzado  : %.4f m\n\n', max(r_cil))

%% 3. GRAFICAS

theta_circ = linspace(0, 2*pi, 300);

% -------------------------------------------------------------------
% Figura 1: Animacion 2D con parpadeo de cables
% -------------------------------------------------------------------
fig1 = figure('Color','w','Name','Fig1 - Animacion hibrida');
hold on; grid on; axis equal
lim = R_cargas * 1.18;
xlim([-lim, lim]); ylim([-lim, lim])
title('Confinamiento hibrido F = q(E + v \times B)  —  RK4', 'FontSize', 11)
xlabel('x (m)'); ylabel('y (m)')

% Circulos de referencia
plot(R_cargas*cos(theta_circ), R_cargas*sin(theta_circ), ...
     'r-',  'LineWidth', 2,   'DisplayName', 'Anillo cargas/cables (R)')
plot(R_limit *cos(theta_circ), R_limit *sin(theta_circ), ...
     'm--', 'LineWidth', 1.5, 'DisplayName', 'R_{limit} (activa rescate)')

% Electrodos / cables (circulo grande, 3 colores segun control)
hElec = gobjects(1, N_cables);
for i = 1:N_cables
    hElec(i) = plot(xC(i), yC(i), 'o', ...
        'MarkerSize', 16, 'LineWidth', 2, ...
        'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'k');
end

% Trayectoria y particula
hTraj = plot(nan, nan, 'b-', 'LineWidth', 0.7, 'DisplayName', 'Trayectoria');
hPart = plot(nan, nan, 'bo', 'MarkerFaceColor','b','MarkerSize', 7);

legend('Anillo (R)', 'R_{limit}', 'Location', 'northeastoutside')

frame_skip = max(1, floor(steps / 400));  % ~400 cuadros en el GIF
primero = true;

for n = 1:frame_skip:steps
    [~, active] = get_currents(traj(n,:), p);
    update_wire_colors(hElec, active);
    set(hTraj, 'XData', traj(1:n,1), 'YData', traj(1:n,2));
    set(hPart, 'XData', traj(n,1),   'YData', traj(n,2));
    drawnow
    if make_gif
        primero = grabar_frame(fig1, gif_file, primero, 0.05);
    end
end
% Dibujar trayectoria completa al final
set(hTraj, 'XData', traj(:,1), 'YData', traj(:,2));
drawnow
if make_gif, fprintf('GIF guardado: %s\n\n', gif_file); end

% -------------------------------------------------------------------
% Figura 2: Radio cilindrico vs tiempo
% -------------------------------------------------------------------
figure('Color','w','Name','Fig2 - Radio vs Tiempo')
hold on; grid on
plot(t*1e9, r_cil, 'b-', 'LineWidth', 1, 'DisplayName','|r(t)|')
yline(R_limit,  'm--', 'LineWidth', 1.5, 'DisplayName', 'R_{limit} (rescate)')
yline(R_escape, 'r--', 'LineWidth', 1.5, 'DisplayName', 'R_{escape} (anillo)')
xlabel('Tiempo (ns)'); ylabel('|r| (m)')
title('Distancia al centro vs tiempo')
legend('Location','best')

% -------------------------------------------------------------------
% Figura 3: Energias vs tiempo
% -------------------------------------------------------------------
figure('Color','w','Name','Fig3 - Energias vs Tiempo')
hold on; grid on
plot(t*1e9, Ek,   'r-', 'LineWidth', 1,   'DisplayName', 'E_k cinetica')
plot(t*1e9, Ep,   'c-', 'LineWidth', 1,   'DisplayName', 'E_p electrica (q·V)')
plot(t*1e9, Etot, 'k-', 'LineWidth', 2.5, 'DisplayName', 'E_{total} = E_k + E_p')
yline(Etot(1), 'k--', 'LineWidth', 1, 'DisplayName', 'E_{total,0} (referencia)')
xlabel('Tiempo (ns)'); ylabel('Energia (J)')
title({'Energias vs tiempo', ...
    'E_{total} conservada (B no hace trabajo) | E_k y E_p intercambian con el campo E'})
legend('Location','best')

% -------------------------------------------------------------------
% Figura 4: Deriva relativa de la energia total (error numerico RK4)
% -------------------------------------------------------------------
figure('Color','w','Name','Fig4 - Validacion numerica')
plot(t*1e9, dE, 'k-', 'LineWidth', 1)
hold on; grid on
yline(0, 'r--', 'LineWidth', 1.5)
xlabel('Tiempo (ns)'); ylabel('\Delta E_{total} / E_{total,0}  (%)')
title('Deriva relativa de energia total — validacion RK4  (ideal \approx 0 %)')

% -------------------------------------------------------------------
% Figura 5: Comparacion trayectoria con B vs sin B
% Demuestra que el campo magnetico SI altera la trayectoria aunque
% no cambie la energia total.
% -------------------------------------------------------------------
fprintf('Corriendo simulacion de referencia sin campo B...\n')
p_sinB        = p;
p_sinB.Bz_0   = 0;      % sin campo uniforme
p_sinB.I_idle   = 0;    % sin corriente en reposo
p_sinB.I_rescue = 0;    % sin corriente de rescate
[traj_sinB, ~, ~, ~] = simular(r0, v0, dt, steps, p_sinB);

figure('Color','w','Name','Fig5 - Comparacion con/sin campo B')

subplot(1,2,1)
hold on; grid on; axis equal
xlim([-R_cargas*1.1, R_cargas*1.1]); ylim([-R_cargas*1.1, R_cargas*1.1])
plot(R_cargas*cos(theta_circ), R_cargas*sin(theta_circ), 'r-',  'LineWidth', 1.5)
plot(R_limit *cos(theta_circ), R_limit *sin(theta_circ), 'm--', 'LineWidth', 1)
plot(xC, yC, 'ko', 'MarkerFaceColor',[0.3 0.3 0.3],'MarkerSize',10)
plot(traj_sinB(:,1), traj_sinB(:,2), 'r-', 'LineWidth', 0.6)
plot(r0(1), r0(2), 'go', 'MarkerFaceColor','g', 'MarkerSize', 8)
title('Solo campo electrico (Bz = 0)')
xlabel('x (m)'); ylabel('y (m)')
legend('Anillo','R_{limit}','Cables','Trayectoria','r_0','Location','northeastoutside')

subplot(1,2,2)
hold on; grid on; axis equal
xlim([-R_cargas*1.1, R_cargas*1.1]); ylim([-R_cargas*1.1, R_cargas*1.1])
plot(R_cargas*cos(theta_circ), R_cargas*sin(theta_circ), 'r-',  'LineWidth', 1.5)
plot(R_limit *cos(theta_circ), R_limit *sin(theta_circ), 'm--', 'LineWidth', 1)
plot(xC, yC, 'ko', 'MarkerFaceColor',[0.3 0.3 0.3],'MarkerSize',10)
plot(traj(:,1), traj(:,2), 'b-', 'LineWidth', 0.6)
plot(r0(1), r0(2), 'go', 'MarkerFaceColor','g', 'MarkerSize', 8)
title(sprintf('Hibrido: E + B  (Bz_0 = %.0e T)', Bz_0))
xlabel('x (m)'); ylabel('y (m)')
legend('Anillo','R_{limit}','Cables','Trayectoria','r_0','Location','northeastoutside')

sgtitle({'Efecto del campo magnetico sobre la trayectoria', ...
    'Mismas condiciones iniciales — la diferencia es solo el campo B'}, ...
    'FontSize', 11)

%% FUNCIONES LOCALES

% ----- Bucle principal de simulacion con RK4 ------------------------------
function [traj, Ek, Ep, t_cpu] = simular(r0, v0, dt, steps, p)
    r = r0; v = v0;
    traj = zeros(steps, 2);  traj(1,:) = r0;
    Ek   = zeros(steps, 1);  Ek(1)     = 0.5 * p.m * sum(v0.^2);
    Ep   = zeros(steps, 1);  Ep(1)     = compute_Ep(r0, p);
    tic
    for n = 1:steps-1
        [r, v]    = rk4_step(r, v, dt, p);
        traj(n+1,:) = r;
        Ek(n+1)   = 0.5 * p.m * sum(v.^2);
        Ep(n+1)   = compute_Ep(r, p);
    end
    t_cpu = toc;
end

% ----- RK4 sobre estado s = [x, y, vx, vy] --------------------------------
function [r_new, v_new] = rk4_step(r, v, dt, p)
    s  = [r, v].';
    k1 = state_deriv(s,               p);
    k2 = state_deriv(s + 0.5*dt*k1,  p);
    k3 = state_deriv(s + 0.5*dt*k2,  p);
    k4 = state_deriv(s +     dt*k3,  p);
    s_new = s + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
    r_new = s_new(1:2).';
    v_new = s_new(3:4).';
end

function ds = state_deriv(s, p)
    r = s(1:2).';
    v = s(3:4).';
    [ax, ay] = compute_accel(r, v, p);
    ds = [v(1); v(2); ax; ay];
end

% ----- Aceleracion total: F = q(E + v x B) / m ----------------------------
function [ax, ay] = compute_accel(r, v, p)

    % == ETAPA 1: campo electrico de Coulomb ==
    % Suma de N cargas puntuales positivas en el anillo.
    % E apunta desde cada carga hacia el punto de campo (repulsion).
    Ex = 0;  Ey = 0;
    for i = 1:p.N_cargas
        dx   = r(1) - p.xQ(i);
        dy   = r(2) - p.yQ(i);
        dist = max(sqrt(dx^2 + dy^2), 1e-9);   % regularizacion
        fac  = p.K * p.Q(i) / dist^3;
        Ex   = Ex + fac * dx;
        Ey   = Ey + fac * dy;
    end

    % == ETAPA 2: campo magnetico uniforme Bz_0 ==
    % Campo constante en direccion z: provoca giros de Larmor en el plano XY.
    Bz = p.Bz_0;

    % == ETAPA 3: campo activo Biot-Savart (control sectorizado) ==
    % Si |r| > R_limit, 3 cables se refuerzan. Su campo se suma a Bz.
    % Convencion: cables pares -> corriente +I, impares -> corriente -I
    % (alternado, igual que en controlDinamicoMag.m).
    I_cables = get_currents(r, p);
    for i = 1:p.N_cables
        dx   = r(1) - p.xC(i);
        dy   = r(2) - p.yC(i);
        dist = max(sqrt(dx^2 + dy^2), 1e-5);
        if mod(i, 2) == 0
            corriente =  I_cables(i);
        else
            corriente = -I_cables(i);
        end
        Bz = Bz + (p.mu0 * corriente) / (2*pi*dist);
    end

    % == Fuerza de Lorentz total ==
    % F = q(E + v x B),  con B = Bz * z_hat
    % v x B = (vy*Bz) x_hat + (-vx*Bz) y_hat
    ax = (p.q / p.m) * (Ex + v(2) * Bz);
    ay = (p.q / p.m) * (Ey - v(1) * Bz);
end

% ----- Energia potencial electrica: Ep = q * V(r) -------------------------
function Ep = compute_Ep(r, p)
    V = 0;
    for i = 1:p.N_cargas
        dx   = r(1) - p.xQ(i);
        dy   = r(2) - p.yQ(i);
        dist = max(sqrt(dx^2 + dy^2), 1e-9);
        V    = V + p.K * p.Q(i) / dist;
    end
    Ep = p.q * V;   % [J]
end

% ----- Control activo: devuelve vector de corrientes I[1..N_cables] --------
function [I, active] = get_currents(r, p)
    I      = ones(1, p.N_cables) * p.I_idle;
    active = struct('on', false, 'idx', [], 'neigh', []);

    if norm(r) > p.R_limit
        ang_p = atan2(r(2), r(1));
        if ang_p < 0, ang_p = ang_p + 2*pi; end

        [~, idx] = min(abs(p.ang_C - ang_p));     % cable mas cercano
        izq      = mod(idx - 2, p.N_cables) + 1;  % vecino izquierdo
        der      = mod(idx,     p.N_cables) + 1;  % vecino derecho

        I(idx) = p.I_rescue;
        I(izq) = p.I_rescue * 0.5;
        I(der) = p.I_rescue * 0.5;

        active.on    = true;
        active.idx   = idx;
        active.neigh = [izq, der];
    end
end

% ----- Actualizar color de los cables en la animacion ---------------------
function update_wire_colors(h, active)
    % Gris  = reposo | Naranja = vecino de rescate | Rojo = cable principal
    set(h, 'MarkerFaceColor', [0.3 0.3 0.3]);
    if active.on
        set(h(active.idx),   'MarkerFaceColor', 'r');
        set(h(active.neigh), 'MarkerFaceColor', [1 0.5 0]);
    end
end

% ----- Guardar cuadro en el GIF ------------------------------------------
function primero = grabar_frame(fig, archivo, primero, retardo)
    try
        fr = getframe(fig);
        [A, map] = rgb2ind(frame2im(fr), 256);
        if primero
            imwrite(A, map, archivo, 'gif', 'LoopCount', Inf, 'DelayTime', retardo);
        else
            imwrite(A, map, archivo, 'gif', 'WriteMode', 'append', 'DelayTime', retardo);
        end
    catch
        exportgraphics(fig, archivo, 'Append', ~primero);
    end
    primero = false;
end
