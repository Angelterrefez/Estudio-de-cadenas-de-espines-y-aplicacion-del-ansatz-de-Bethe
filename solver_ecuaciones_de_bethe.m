% -*- coding: utf-8 -*-
% -*- babel: spanish -*-

% -------------------------------------------------------------------------
% Creado el 03/07/2025
% @author: José Ángel Terrero Fernández
% -------------------------------------------------------------------------

%% Bethe ansatz "solver" para la cadena XXX espín-1/2 (Periódica, B=0)
% -------------------------------------------------------------------------
% El script resuelve numéricamente las ecuaciones de Bethe a partir del 
% estado introducido en función de N, M y los m_j. Obtiene los parámetros 
% λ_j, ajusta los cuasimomentos k_j y calcula la energía del estado.
% Este script está optimizado para soluciones con dos magnones (m_1,m_2)
% -------------------------------------------------------------------------
% Se pide establecer la configuración del sistema (N,M) y el 
% estado seleccionado (los m_j tomando valores 0,...,N-1).
% La rutina que sigue el script es en primer lugar resolver las ecuaciones 
% de Bethe en λ_j y luego ajustar la rama logarítmica para dar los k_j. Por
% último determina la energía del sistema. A continuación se presentan las
% ecuaciones empleadas por el código.
%
% Ecuaciones de Bethe en función de las rapidities (λ_j):
%   ((λ_j + i)/(λ_j - i))^N = ∏_{l≠j} ((λ_j - λ_l + 2i)/(λ_j - λ_l - 2i)),
%   j = 1,...,M.
%
% Conversión de λ_j a k_j ajustando la rama:
%   k_j0 = 2*atan(1/λ_j),
%   k_j = k_j0 + 2π * round((2π*m_j/N - k_j0)/(2π)).
%
% Energía:
%   E = ∑_{j=1}^M 2/(λ_j^2 + 1).
% -------------------------------------------------------------------------

function bethe_solver()
    % Limpia la ventana de comandos y muestra el encabezado.
    clc;
    fprintf('===============================\n');
    fprintf(' Solucionador Bethe Ansatz XXX  \n');
    fprintf('===============================\n');
    
    % =====================================================================
    % Entrada de parámetros del sistema
    % =====================================================================
    % N: Número de espines en la cadena (debe ser entero positivo).
    % M: Número de magnones/espines volteados.
    N = input('Introduce el número de espines (N): ');
    M = input('Introduce el número de magnones (M): ');
    
    % =====================================================================
    % Lectura de números cuánticos m_j
    % =====================================================================
    % Los m_j definen el estado cuántico y deben ser enteros en [0, N-1].
    fprintf('\nIntroduce los m_j que definen al estado (0,...,N-1):\n');
    m = zeros(M,1);
    for j = 1:M
        m(j) = input(sprintf('m_(%d) = ', j));
    end
    
    % Validación de los valores de m_j.
    if any(m<0) || any(m> N-1) || any(floor(m)~=m)
        error('Los m_j deben ser enteros en [0,N-1]');
    end

    % =====================================================================
    % Configuración inicial para resolver ecuaciones de Bethe
    % =====================================================================
    % Función anónima que define las ecuaciones de Bethe (F(λ) = 0).
    fun = @(lam) bethe_eqs(lam, N, M);
    
    % Semilla inicial para λ_j (aproximación para sistema no
    % interactuante).
    % Fórmula: λ_j ≈ 0.5 * cot(k_j/2). Relación del límite termodinámico.
    lambda0 = 0.5 * cot(pi * m / N);
    
    % Manejo de singularidades: Si m_j = 0 o N, cot(π*m_j/N) es infinito.
    bad = ~isfinite(lambda0);  % Detecta valores NaN/Inf.
    if any(bad)
        warning('Semilla singular, usando valores aleatorios.');
        lambda0(bad) = 1e6;  % Valor grande para λ, que se traduce en un
                             % infinito aproximado.
    end

    % =====================================================================
    % Resolución numérica con fsolve
    % =====================================================================
    % Opciones del solver: Tolerancias ajustadas para precisión.
    opts = optimoptions('fsolve', ...
        'Display', 'off', ...      % Sin salida iterativa
        'TolFun', 1e-12, ...       % Tolerancia en función objetivo
        'TolX', 1e-12, ...         % Tolerancia en variables
        'MaxFunEvals', 1e4, ...    % Máximas evaluaciones de función
        'MaxIter', 1e4);           % Máximas iteraciones
    
    % Llamada al solver de ecuaciones no lineales.
    [lambda_sol, ~, exitflag] = fsolve(fun, lambda0, opts);
    
    % Verifica convergencia del solver.
    if exitflag <= 0
        warning('No convergió fsolve en λ (exitflag=%d).', exitflag);
    end
    
    % =====================================================================
    % Resultados: Rapidities λ_j
    % =====================================================================
    fprintf('\n-------------------------------------');
    fprintf('\nRapidities (λ_j):\n');
    disp(lambda_sol(:));  % Muestra los λ_j como vector columna.

    % =====================================================================
    % Cálculo de cuasimomentos k_j con corrección de rama
    % =====================================================================
    % k_j0: Valor inicial sin ajuste de rama (puede estar en rama 
    % equivocada).
    k0 = 2 * atan(1 ./ lambda_sol);  % Relación λ_j = 1/(tan(k_j/2))
    
    % k_target: Valores esperados de k_j según números cuánticos m_j.
    k_target = 2*pi*m / N;  % Cuantización de momento en cadena periódica.
    
    % Ajuste de rama: Añade múltiplos de 2π para alinear k_j con k_target.
    k = k0 + 2*pi * round((k_target - k0) / (2*pi));
    
    fprintf('\n-------------------------------------');
    fprintf('\nCuasimomentos ajustados (k_j):\n');
    disp(k(:));  % Muestra los k_j corregidos.

    % =====================================================================
    % Cálculo de la energía del estado
    % =====================================================================
    % Energía por magnón: 2/(λ_j^2 + 1). Fórmula del modelo XXX.
    E = sum(2 ./ (lambda_sol.^2 + 1));
    fprintf('\n-------------------------------------');
    fprintf('\nEnergía E = %g\n', E);
    fprintf('-------------------------------------\n');
end

% =========================================================================
% Función: bethe_eqs
% =========================================================================
% Define las ecuaciones de Bethe para el modelo XXX.
% Entradas:
%   lambda: Vector de rapidities λ_j
%   N: Número de espines
%   M: Número de magnones
% Salida:
%   F: Vector donde F(j) = 0 es la ecuación a resolver en λ_j.
% =========================================================================
function F = bethe_eqs(lambda, N, M)
    i = 1i;  % Unidad imaginaria.
    F = zeros(M,1);  % Inicializa vector de ecuaciones.
    
    % Bucle sobre cada magnón j.
    for j = 1:M
        % Lado izquierdo (LHS): Término de periodicidad.
        lhs = ((lambda(j) + i) / (lambda(j) - i))^N;
        
        % Lado derecho (RHS): Productorio de interacciones magnón-magnón.
        rhs = 1;
        for l = 1:M
            if l ~= j
                rhs = rhs * ((lambda(j) - lambda(l) + 2*i) / ...
                             (lambda(j) - lambda(l) - 2*i));
            end
        end
        
        % Ecuación j-ésima: LHS - RHS = 0.
        F(j) = lhs - rhs;
    end
end

% Ejecuta la función principal.
bethe_solver();
