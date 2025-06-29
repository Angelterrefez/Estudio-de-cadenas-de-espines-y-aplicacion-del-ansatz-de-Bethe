% -*- coding: utf-8 -*-
% -*- babel: spanish -*-

% -------------------------------------------------------------------------
% Creado el 3/07/2025
% @author: José Ángel Terrero Fernández
% -------------------------------------------------------------------------

%% Código para comparar la solución proporcionada por el ansatz de Bethe 
%% con la diagonalización numérica del hamiltoniano.
% -------------------------------------------------------------------------
% Este script construye el hamiltoniano del modelo de Heisenberg XXX 
% para N espines con condiciones de contorno periódicas, incluyendo un 
% campo magnético externo B; y lo diagonaliza, calculando los autovalores y
% autovectores.
%
% Luego, resuelve las ecuaciones de Bethe en cada subespacio de magnones y 
% calcula la energía propia correspondiente a cada estado hasta replicar 
% todas las soluciones numéricas resultado de la diagonalización. Se 
% imprime la combinación de números cuánticos de Bethe (m_j) y los 
% cuasimomentos de los magnones que definen al autoestado.
%
% Finalmente, representa ambas soluciones de autovalores en una misma
% gráfica de energía frente a M (número de espines volteados).
% -------------------------------------------------------------------------
clear; close all; clc;

disp(['Modelo cadena de espines de Heisenberg XXX con' ...
    ' condiciones de contorno periódicas y campo externo aplicado (B).']);
disp(' ');

% Parámetros del sistema
N = input('Introduzca el número de espines (N): ');
J = input('Introduzca el valor de la constante de intercambio (J): ');
B = input('Introduzca la intensidad del campo magnético (B): ');
hbar = 1; % hbar: Constante reducida de Planck

%% 1. Diagonalización numérica del hamiltoniano. 
disp(['-----------------------------------------------------------------' ...
    '-------------------------']);
disp(['1. Cálculo de autovalores y autovectores del hamiltoniano mediante ' ...
    'diagonalización numérica.']);
disp(['-----------------------------------------------------------------' ...
    '-------------------------']);
[eigenvalues_diag, M_diag, eigenvectors] = diagonalize_XXX(N, J, B);

% Mostrar resultados
disp('Autovalores:');
disp(eigenvalues_diag');
disp('Matriz de autovectores:');
disp(eigenvectors);

%% 2. Calcular energías con ansatz de Bethe para todos los M.
disp(['-----------------------------------------------------------------' ...
    '-------------------------']);
disp('2. Soluciones dadas por el ansatz de Bethe.');
disp(['-----------------------------------------------------------------' ...
    '-------------------------']);

E_bethe = cell(N+1, 1);  % Almacenar energías por sector M
E0 = -J*N/4;  % Energía del estado de referencia

for M_val = 0:N

    fprintf('\n--------------------------------------------------\n');    
    fprintf('Para M = %d hay %d estados propios\n', M_val, sum(M_diag == M_val));
    fprintf('--------------------------------------------------\n');

    % Caso M=0 (estado de referencia)
    if M_val == 0
        E_bethe{1} = E0 - B*(N/2 - 0);
        
        fprintf('\n--------------------------------------------------\n');
        fprintf('Estado para M = %d\n', M_val);
        fprintf('Energía E = %.12f\n', E_bethe{1});
        fprintf('--------------------------------------------------\n');
        continue;
    end
    
    % Caso M=N (estado todos los espines volteados)
    if M_val == N
        E_bethe{N+1} = E0 - B*(N/2 - N);

        fprintf('\n--------------------------------------------------\n');
        fprintf('Estado para M = N = %d\n', M_val);
        fprintf('Energía E = %.12f\n', E_bethe{N+1});
        fprintf('--------------------------------------------------\n');        
        continue;
    end
    
    % Caso general M_val (1 a N-1)
    if M_val == 1
        % Solución analítica para 1 magnón
        k_values = 2*pi*(0:N-1)/N;  % Posibles momentos
        energies = zeros(1, N);
        
        for m1 = 0:N-1
            k1 = 2*pi*m1/N;
            % Energía para este magnón
            en = E0 + J*(1 - cos(k1));
            energies(m1+1) = en;
            fprintf('\n--------------------------------------------------\n');
            fprintf('Estado para M = %d\n', M_val);
            fprintf('Número cuántico m_1 = [ %d ]\n', m1);
            fprintf('Cuasimomento k_1 = [ %s ]\n', num2str(k1, '%.6f '));
            fprintf('Energía E = %.12f\n', en);
            fprintf('--------------------------------------------------\n');
        end
        
        % Eliminar duplicados en el círculo de momentos
        energies = unique(round(energies*1e8)/1e8); % Elimina duplicados numéricos
        E_bethe{M_val+1} = energies - B*(N/2 - M_val);

    else
        % Para M_val > 1, usar el solver de Bethe para cada combinación de m_j
        energies = [];
        tol_match = 1e-1;  % Tolerancia para coincidencia de energía via el ansatz con resultado de diagonalización
        eigs_diag_M = real(eigenvalues_diag(M_diag == M_val));  % Autovalores para este M_val según diagonalización

        % Generar todas las combinaciones posibles de m_j
        m_vals = 0:N-1;
        grid = cell(1, M_val);
        [grid{:}] = ndgrid(m_vals);

        % Cada fila es una combinación
        m_combinations = reshape(cat(M_val+1, grid{:}), [], M_val);
        
        % Sin importar el orden
        m_combinations = unique(sort(m_combinations, 2), 'rows');

        soluciones_para_M = [];

        % Bucle que resuelve las ecuaciones para cada posible combinación de los m_j 
        for idx = 1:size(m_combinations, 1)
            m_j = m_combinations(idx, :);
            match_found = false;
            attempts = 0;
    
            % Se resuelven las ecuaciones de Bethe hasta encontrar el
            % autoestado (definido por la diagonalización). Se para si a
            % las 100 iteraciones no se ha logrado.
            while ~match_found && attempts < 100
                attempts = attempts + 1;

            % Condición inicial para k_j, evitando singularidades
            epsilon = 1e-3;
            k0 = mod(2*pi*(m_j)/N + epsilon * randn(1, M_val), 2*pi); % Perturbación aleatoria pequeña

            % Verificación previa
            try
                test_output = bethe_eqs_full(k0, m_j, N);
                if any(isnan(test_output)) || any(isinf(test_output))
                    continue;  % Saltar combinaciones malas
                end
            catch
                continue;  % Saltar si ocurre error en evaluación
            end
            % Resolver las ecuaciones de Bethe para la combinación de m_j
            options = optimoptions('fsolve', ...
                'Display', 'off', ...
                'TolFun', 1e-12, ...
                'TolX', 1e-12, ...
                'MaxFunEvals', 1e4, ...
                'MaxIter', 1e4);
            [k_sol, ~, exitflag] = fsolve(@(k) bethe_eqs_full(k, m_j, N), k0, options);

            if exitflag <= 0 || any(isnan(k_sol)) || any(isinf(k_sol))
                continue;
            end

            % Ordenar y evitar duplicados
            k_sol = sort(mod(k_sol, 2*pi));
            if ~isempty(soluciones_para_M)
                if any(all(abs(soluciones_para_M - k_sol) < 1e-6, 2))
                    continue;
                end
            end
            soluciones_para_M = [soluciones_para_M; k_sol];

            % Calcular energía
            EM = E0 + J * hbar^2 * sum(1 - cos(k_sol)) - B * (N/2 - M_val) * hbar;

            % Comparar con autovalores numéricos
            if any(abs(EM - eigs_diag_M) < tol_match)
                match_found = true;

                energies = [energies, EM];


% Mostrar resultados
fprintf('\n--------------------------------------------------\n');
fprintf('Estado para M = %d\n', M_val);
fprintf('Combinación m_j = [ %s ]\n', num2str(m_j));
fprintf('Cuasimomentos k_j = [ %s ]\n', num2str(k_sol, '%.6f '));
fprintf('Energía E = %.12f\n', EM);
fprintf('--------------------------------------------------\n');
            end
            end
            
        end

        % Eliminar duplicados numéricos
        energies = unique(round(energies*1e8)/1e8);
        E_bethe{M_val+1} = energies - B*(N/2 - M_val);
    end
end

%% 3. Graficar resultados
figure; hold on;

for M_val = 0:N
    % Diagonalización (punto azul)
    mask = (M_diag == M_val);
    h_diag = scatter(M_val*ones(sum(mask),1), real(eigenvalues_diag(mask)), ...
        160,[0.2 0.6 1], 'filled', 'MarkerEdgeColor', 'k', 'SizeData', 160);
    
    % Bethe (equis roja)
    if ~isempty(E_bethe{M_val+1})
        h_bethe = scatter(M_val*ones(numel(E_bethe{M_val+1}),1), ...
            real(E_bethe{M_val+1}),180,'rx','LineWidth',1.5,'SizeData', 300);
    end
end

legend([h_diag, h_bethe], ...
       {'Diagonalización', 'Bethe ansatz'}, ...
       'Location', 'best', 'FontSize', 18);

% Configuración del gráfico
xlabel('Número de magnones, M', 'FontSize', 22);
ylabel('Energía, E_M', 'FontSize', 22);
title(sprintf(['Comparación diagonalización y soluciones de Bethe ' ...
    '(N=%d, B=%d)'], N, B), 'FontSize', 24);
xticks(0:1:N);
set(gca,'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 20);
xlim([min(M_diag), max(M_diag)]);
ylim([min(eigenvalues_diag)-0.5, max(eigenvalues_diag) + 0.5]);
hold off;

%% ------------------------------------------------------------------------
% Función: Diagonalización del hamiltoniano XXX.
function [eigenvalues, M_vector, V] = diagonalize_XXX(N, J, B)
    % Operadores de Pauli
    sx = [0 1; 1 0]/2;
    sy = [0 -1i; 1i 0]/2;
    sz = [1 0; 0 -1]/2;
    id = eye(2);

    % Inicializar la matriz hamiltoniana (tamaño 2^N x 2^N, ya que cada 
    % espín tiene 2 configuraciones posibles)
    dim = 2^N;
    H = zeros(dim);

% =====================================================================
% Construcción del hamiltoniano
% =====================================================================
% El hamiltoniano incluye:
% - Términos de interacción entre primeros vecinos:
% (-J(sx_j sx_{j+1} + sy_j sy_{j+1} + sz_j sz_{j+1}))
% - Condiciones de contorno periódicas: s_{j+N} = s_j
% - Término del campo magnético externo: (-B Σ sz_j)
% ---------------------------------------------------------------------
    % Construir términos de interacción
    for i = 1:N-1
        for op = {sx, sy, sz}
            term = 1;
            for j = 1:N
                if j == i || j == i+1
                    term = kron(term, op{1});
                else
                    term = kron(term, id);
                end
            end
            H = H - J*term;
        end
    end
    
    % Término condiciones de contorno periódicas
    if N > 1
        for op = {sx, sy, sz}
            term = 1;
            for j = 1:N
                if j == 1 || j == N
                    term = kron(term, op{1});
                else
                    term = kron(term, id);
                end
            end
            H = H - J*term;
        end
    end
    
    % Término de campo magnético
    for i = 1:N
        term = 1;
        for j = 1:N
            if j == i
                term = kron(term, sz);
            else
                term = kron(term, id);
            end
        end
        H = H - B*term;
    end
    
    % Diagonalización
    [V, D] = eig(H);
    eigenvalues = diag(D);
    
    % Calcular M (espines volteados)
    M_vector = calculate_M_from_Sz(V, N);
end

%% ------------------------------------------------------------------------
% Función: Diagonaliza S_T^z en la misma base para cada autovector y le
% asigna el valor de M correspondiente.
function M_vector = calculate_M_from_Sz(V, N)
    % Calcula M como el número de espines volteados usando S_z^total
    sz = [1 0; 0 -1]/2;
    id = eye(2);
    Sz_total = zeros(2^N);

    % Construir el operador S^z total
    for i = 1:N
        term = 1;
        for j = 1:N
            if j == i
                term = kron(term, sz);
            else
                term = kron(term, id);
            end
        end
        Sz_total = Sz_total + term;
    end

    % Para cada autovector, calcula el valor esperado de S^z total
    M_vector = zeros(size(V, 2), 1);  % M para cada autovector
    for i = 1:size(V, 2)
        psi = V(:, i);
        sz_total_exp = real(psi' * Sz_total * psi);
        % M = N/2 - ⟨S^z⟩
        M_vector(i) = round(N/2 - sz_total_exp);  % redondear a entero
    end
end

%% ------------------------------------------------------------------------
% Función: Resolver ecuaciones de Bethe para una combinación de m_j
function F = bethe_eqs_full(k, m, N)
    M = length(k);
    F = zeros(1, M);
    for j = 1:M
        sum_theta = 0;
        for l = 1:M
            if l ~= j
                delta_cot = cot_safe(k(j)/2) - cot_safe(k(l)/2);
                theta_jl = 2 * atan2(1, 0.5 * delta_cot);
                sum_theta = sum_theta + theta_jl;
            end
        end
        F(j) = N * k(j) - 2 * pi * m(j) - sum_theta;
    end
end

%% ------------------------------------------------------------------------
% Función: Evita singularidades
function c = cot_safe(x)
    tol = 1e-6;
    xmod = mod(x, pi);
    if abs(xmod) < tol || abs(xmod - pi) < tol
        c = 1 / tan(x + tol * sign(x + 1e-12)); 
    else
        c = cot(x);
    end
end
