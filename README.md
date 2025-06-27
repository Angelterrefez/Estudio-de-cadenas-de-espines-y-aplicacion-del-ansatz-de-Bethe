# Estudio de cadenas de espines y aplicación del ansatz de Bethe.

Autor: José Ángel Terrero Fernández.

Códigos para replicar los resultados numéricos obtenidos en un TFG del Grado en Física.
Cadena unidimensional de espines de Heisenberg, caso isótropo XXX con condiciones de contorno periódicas a campo magnético externo uniforme B. Caso ferromagnético.

Departamento de Física Atómica, Molecular y Nuclear. Facultad de Física. Universidad de Sevilla.
Curso Académico 2024/25.

# Archivos .m
1. espectro_energetico.m: grafica el espectro energético obtenido con el ansatz de Bethe para N espines y M magnones.
2. solver_ecuaciones_de_bethe.m: resuelve las ecuaciones de Bethe en su correspondiente subespacio de magnones a partir del estado seleccionado con los números cuánticos para N espines y M magnones. Determina los cuasimomentos de los magnones y la energía del estado.
3. diagonalizacion_vs_bethe.m: compara las soluciones analíticas obtenidas via el ansatz de Bethe con las de la diagonalización numérica del hamiltoniano para N espines y M magnones.
