% -*- coding: utf-8 -*-
% -*- babel: spanish -*-

% -------------------------------------------------------------------------
% Creado el 03/07/2025
% @author: José Ángel Terrero Fernández
% -------------------------------------------------------------------------

%% Script para representar el continuo de M-magnones, el estado ligado 
%% y la dispersión de un solo magnón con campo externo magnético aplicado.
% -------------------------------------------------------------------------
% Este script calcula y grafica el espectro de energía de una cadena de 
% espines bajo un campo magnético externo B, J=hbar=1, incluyendo:
%   1. Dispersión de un solo magnón.
%   2. Estado ligado de M magnones (solución tipo cuerda de Bethe).
%   3. Continuo de energía en límite termodinámico.
%   4. Soluciones de las ecuaciones de Bethe para M=2.
% Estados discretos para cadenas finitas para M>2. 
% -------------------------------------------------------------------------

clear; close all; clc;

%% Parámetros del sistema
N = input('Número de espines (N): ');
M = input('Número de magnones (M): ');
B = input('Intensidad del campo magnético aplicado (B): ');

%% 1. Dispersión de un solo magnón (M>=1)
% -------------------------------------------------------------------------
% Energía de un solo magnón (M=1): E_1 = 1 - cos(k) - B*(N/2 - 1).
% -------------------------------------------------------------------------
Kvec = linspace(-pi, pi, 100); % K ∈ [-π, π]

if M >= 1
    E_1 = 1 - cos(Kvec) - B*(N/2 - 1);
end

%% 2. Estado ligado (M>=2)
% -------------------------------------------------------------------------
% Los estados ligados tienen energía por debajo del continuo.
% Para el modelo XXX, la energía del estado ligado es: 
% E_bound = (1/M)(1 - cos(K)) - B*(N/2 - M).
% -------------------------------------------------------------------------
if M >= 2
    E_bound = (1/M)*(1 - cos(Kvec)) - B*(N/2 - M);
end

%% 3. Continuo de energía para M magnones
% -------------------------------------------------------------------------
% Para M=2, en el límite termodinámico, el continuo va desde:
% E_min = 2(1 - cos(K/2))- B*(N/2 - M), a 
% E_max = 2(1 + cos(K/2))- B*(N/2 - M).
% -------------------------------------------------------------------------
if M==2
    E_min_th = 2 * (1 - cos(Kvec/2)) -B*(N/2 - M);
    E_max_th = 2 * (1 + cos(Kvec/2)) -B*(N/2 - M);

% -------------------------------------------------------------------------
% En el límite termodinámico, el espectro forma un continuo entre energías 
% mínima y máxima, pero no hay una expresión explícita. Se muestrean 
% configuraciones aleatorias de momentos para estimar el continuo si M>2.
% -------------------------------------------------------------------------    
elseif M > 2
    fprintf('Calculando continuo con muestreo aleatorio para M = %d...\n',M);
    N_samples = 5000; % Número de configuraciones aleatorias por K

    E_M_min = zeros(size(Kvec));
    E_M_max = zeros(size(Kvec));
    
    for idx = 1:length(Kvec)
        K = Kvec(idx);
        % Genera M-1 momentos aleatorios (distribución uniforme en [-π, π]).
        rand_momenta = -pi + 2*pi*rand(N_samples, M-1);
         % Calcula el M-ésimo momento para conservar K = sum(k_i).
        k_last = K - sum(rand_momenta, 2);
        % Ajusta k_last al rango [-π, π] (periodicidad).
        k_last = mod(k_last + pi, 2*pi) - pi;
        % Energía total de la configuración (Hamiltoniano XXX con campo B):
        % E = sum_{j=1}^M (1 - cos(k_j)) - B*(N/2 - M).
        E_tot = sum(1 - cos(rand_momenta), 2) + ...
            (1 - cos(k_last)) - B*(N/2 - M);
        % Registra mínimos y máximos para el continuo.
        E_M_min(idx) = min(E_tot);
        E_M_max(idx) = max(E_tot);
    end
end

%% 4. Estados discretos para cadena finita.
% -------------------------------------------------------------------------
% Para cadenas finitas con N espines y M=2, se resuelven las
% ecuaciones de Bethe.
% -------------------------------------------------------------------------
if M == 2
    fprintf('Resolviendo ecuaciones de Bethe para M=2...\n');
    m = 0:(N-1);
    [m1, m2] = meshgrid(m, m);
    idx = m1 < m2;
    m1 = m1(idx);
    m2 = m2(idx);
    total_pairs = length(m1);
    
    K_scatt = zeros(total_pairs, 1);
    E_scatt = zeros(total_pairs, 1);
    valid_solution = false(total_pairs, 1);
    
    max_iter = 200;
    tol = 1e-10;
    
    for p = 1:total_pairs
        % Inicializar momentos
        k0 = [2*pi*m1(p)/N; 2*pi*m2(p)/N];
        
        % Evitar singularidades
        if any(k0 == 0)
            k0(k0 == 0) = 1e-8;
        end
        
        % Función para calcular las soluciones a las ecuaciones de Bethe.
        bethe_residuals = @(k) [
            N*k(1) - 2*pi*m1(p) - 2*acot(0.5*(cot(k(1)/2) - cot(k(2)/2)));
            N*k(2) - 2*pi*m2(p) + 2*acot(0.5*(cot(k(1)/2) - cot(k(2)/2)))
        ];
        
        options = optimset('Display','off', 'TolFun', tol, ...
            'MaxIter', max_iter);
        try
            [ksol, ~, exitflag] = fsolve(bethe_residuals, k0, options);
            if exitflag > 0 && all(isreal(ksol))
                K_total = sum(ksol);
                K_total = mod(K_total + pi, 2*pi) - pi;
                E_total = sum(1 - cos(ksol)) - B*(N/2 - 2);
                
                K_scatt(p) = K_total;
                E_scatt(p) = E_total;
                valid_solution(p) = true;
            end
        catch
            continue;
        end
    end
    
    % Filtrar soluciones válidas
    K_scatt = K_scatt(valid_solution);
    E_scatt = E_scatt(valid_solution);

% -------------------------------------------------------------------------
% Para cadenas finitas con N espines, si M>2, los momentos permitidos son:
% k_n = 2πn/N, donde n = -floor(N/2), ..., ceil(N/2)-1 (periodicidad).
% -------------------------------------------------------------------------
elseif M >= 2
    n_vals = -floor(N/2):(ceil(N/2)-1); % Índices cuánticos
    k_allowed = 2*pi*n_vals/N; % Momentos permitidos
    num_allowed = length(k_allowed);
    combos = nchoosek(1:num_allowed, M);  % Todas las combinaciones de 
                                          % M momentos
    ncomb = size(combos, 1);
    % Manejo computacional: Si hay demasiadas combinaciones, se muestrean:
    max_combos = 1e7; % Umbral para muestreo aleatorio
    if ncomb > max_combos
        fprintf(['El número de combinaciones (%d) es muy alto. ' ...
        'Se muestrearán aleatoriamente 1000000 estados.\n'], ncomb);
        num_samples_finite = min(1000000, ncomb);
        idx_rand = randperm(ncomb, num_samples_finite);
        combos = combos(idx_rand, :);
        ncomb = num_samples_finite;
    end
    % Calcular momentos totales (K) y energías (E) para cada combinación.
    K_finite = zeros(ncomb, 1);
    E_finite = zeros(ncomb, 1);
    for i = 1:ncomb
        k_set = k_allowed(combos(i, :)); % Seleccionar M momentos
        K_tot = sum(k_set); % Momento total
        % Ajustar K_tot al rango [-π, π] (periodicidad)
        K_tot = mod(K_tot + pi, 2*pi) - pi;
        % Calcular energía de la configuración
        E_tot = sum(1 - cos(k_set)) - B*(N/2 - M);
        K_finite(i) = K_tot;
        E_finite(i) = E_tot;
    end
end

%% Graficar resultados
figure; hold on;

% Dispersión de un magnón (línea negra discontinua)
if M >= 1
    plot(Kvec, E_1, 'k--', 'LineWidth', 2);
end

% Estado ligado (curva roja)
if M >= 2
    plot(Kvec, E_bound, 'r-', 'LineWidth', 2);
end

% Continuo de energía (área sombreada en azul)
if M == 2
    fill([Kvec, fliplr(Kvec)], [E_min_th, E_max_th], 'b', ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none');
elseif M > 2
    fill([Kvec, fliplr(Kvec)], [E_M_min, fliplr(E_M_max)], 'b', ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% Estados discretos para cadena finita (puntos azules)
if M==2
    plot(K_scatt,E_scatt,'bo','MarkerFaceColor','b','MarkerSize',4);
elseif M >= 2
    plot(K_finite,E_finite,'bo','MarkerFaceColor','b', 'MarkerSize',4);
end

% Configuración del gráfico
xlabel('Número de onda, k', 'FontSize', 22);
ylabel('Energía, E_{M}-E_{0}', 'FontSize', 22);
title(sprintf(['Espectro de energía para %d espines con %d magnones ' ...
    'a campo externo B = %d'], N, M, B), 'FontSize', 24);
grid off;
xlim([-pi, pi]);
xticks(-pi:pi/6:pi); 
xticklabels({'-\pi', '-5\pi/6', '-2\pi/3', '-\pi/2', '-\pi/3', '-\pi/6', ...
    '0', '\pi/6', '\pi/3', '\pi/2', '2\pi/3', '5\pi/6', '\pi'}); 
set(gca, 'FontSize', 20);
hold off;
