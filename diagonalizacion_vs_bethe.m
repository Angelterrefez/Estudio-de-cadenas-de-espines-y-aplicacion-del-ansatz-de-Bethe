%% Código para comparar la solución proporcionada por el ansatz de Bethe 
%% con la diagonalización numérica del hamiltoniano.
% -------------------------------------------------------------------------
% Este script construye el hamiltoniano del modelo de Heisenberg XXX 
% para N espines con condiciones de contorno periódicas, incluyendo un 
% campo magnético externo B; y lo diagonaliza, calculando los autovalores y
% autovectores.
% Luego, resuelve las ecuaciones de Bethe proporcionadas por el ansatz y
% calcula las energías propias correspondientes a cada estado.
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
    % Caso M=0 (estado de referencia)
    if M_val == 0
        E_bethe{1} = E0 - B*(N/2 - 0);
        % Mostrar resultados
fprintf('\n--------------------------------------------------\n');
fprintf('Estado para M = %d\n', M_val);
fprintf('Energía E = %.12f\n', E_bethe{1});
fprintf('--------------------------------------------------\n');
        continue;
    end
    
    % Caso M=N (estado todos los espines hacia abajo)
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
        % Generar todas las combinaciones posibles de m_j (sin repetición)
        m_combinations = nchoosek(0:N-1, M_val);
        for idx = 1:size(m_combinations, 1)
            m_j = m_combinations(idx, :);
            % Resolver las ecuaciones de Bethe para esta combinación de m_j
            lambda_sol = solve_bethe_eqs(N, M_val, m_j);
            if ~isempty(lambda_sol)
% -------------------------------------------------------------------------
% Conversión de λ_j a k_j ajustando la rama:
%   k_j0 = 2*atan(1/λ_j),
%   k_j = k_j0 + 2π * round((2π*m_j/N - k_j0)/(2π)).
%
% Energía:
%   E = ∑_{j=1}^M 2/(λ_j^2 + 1).
% -------------------------------------------------------------------------
                % Calcular cuasimomentos k_j
                k_j = 2 * atan(1 ./ lambda_sol);
                % Asegurar valores en el rango [0, 2π)
                k_j = mod(real(k_j), 2*pi);

                % Calcular la energía
                E = sum(2 ./ (lambda_sol.^2 + 1));
                En = E0 + J * E;
                energies = [energies, En];

% Mostrar resultados
fprintf('\n--------------------------------------------------\n');
fprintf('Estado para M = %d\n', M_val);
fprintf('Combinación m_j = [ %s ]\n', num2str(m_j));
fprintf('Rapidities λ_j  = [ %s ]\n', num2str(lambda_sol, '%.6f '));
fprintf('Cuasimomentos k_j = [ %s ]\n', num2str(k_j, '%.6f '));
fprintf('Energía E = %.12f\n', En);
fprintf('--------------------------------------------------\n');

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
    M_vector = zeros(size(eigenvalues));
    for i = 1:length(eigenvalues)
        [~, idx] = max(abs(V(:, i)));  % Componente dominante
        state_index = idx - 1;         % Índice 0-base
        bin_state = dec2bin(state_index, N) - '0';  % Convertir a binario
        M_vector(i) = sum(bin_state);  % Contar espines volteados
    end
end

%% ------------------------------------------------------------------------
% Función: Resolver ecuaciones de Bethe para una combinación de m_j
function lambda_sol = solve_bethe_eqs(N, M, m_j)

% =========================================================================
% Resolución de las ecuaciones de Bethe según el estado seleccionado 
% (combinación de m_j)
% =========================================================================
    % Función anónima que define las ecuaciones de Bethe
    fun = @(lam) bethe_eqs(lam, N, M);
    
    % Semilla inicial para λ_j (aproximación para sistema no
    % interactuante).
    % Fórmula: λ_j ≈ M * cot(k_j/2). Relación del límite termodinámico.
    lambda0 = M * cot(pi * m_j / N);
    
    % Manejo de singularidades: Si m_j = 0 o N, cot(π*m_j/N) es infinito.
    bad = ~isfinite(lambda0);  % Detecta valores NaN/Inf.
    if any(bad)
        warning('Semilla singular, usando valores aleatorios.');
        lambda0(bad) = 1e6;  % Valor grande para λ, que se traduce en un
                             % infinito aproximado.
    end
    
    % Opciones del solver
    opts = optimoptions('fsolve', ...
        'Display', 'off', ...
        'TolFun', 1e-12, ...
        'TolX', 1e-12, ...
        'MaxFunEvals', 1e4, ...
        'MaxIter', 1e4);
    
    % Resolver ecuaciones de Bethe
    [lambda_sol, fval, exitflag] = fsolve(fun, lambda0, opts);
    
    % Verificar convergencia
    if exitflag <= 0
        warning(['fsolve no convergió para m_j = [%s], exitflag = %d, ' ...
            '||fval|| = %.2e'], num2str(m_j), exitflag, norm(fval));
        lambda_sol = [];
    end
end

%% ------------------------------------------------------------------------
% Función: Ecuaciones de Bethe
function F = bethe_eqs(lambda, N, M)
% =====================================================================
% Construcción de las ecuaciones de Bethe
% =====================================================================
% Ecuaciones de Bethe en función de las rapidities (λ_j):
%   ((λ_j + i)/(λ_j - i))^N = ∏_{l≠j} ((λ_j - λ_l + 2i)/(λ_j - λ_l - 2i)),
%   j = 1,...,M.
% -------------------------------------------------------------------------
    i = 1i;
    F = zeros(M,1);
    
    for j = 1:M
        lhs = ((lambda(j) + i) / (lambda(j) - i))^N;
        rhs = 1;
        for l = 1:M
            if l ~= j
                rhs = rhs * ((lambda(j) - lambda(l) + 2*i) / ...
                             (lambda(j) - lambda(l) - 2*i));
            end
        end
        F(j) = lhs - rhs;
    end
end