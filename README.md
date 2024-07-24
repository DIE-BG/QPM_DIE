# QPM_DIE

Herramienta de Corrimientos: QPM- Versión 2

Modelo: QPM.

Versión 2: versión DIE con ajustes en la estructura.

Estructura de la herramienta de ejecución:

PreProcessing:
Carga de datos primitivos (en el nivel y frecuencia de la fuente).
Cálculos y transformaciones de frecuencia para gráficas de preprocesamiento y generación de base de datos con variables observables (ver .mod).

Filtrado:
Con base en la base de datos de observables generadas en PreProcessing se utiliza el filtro de Kalman para recuperar el resto de variables y generar la base de datos de condiciones iniciales para la simulación.

Simulación:
Utiliza la base de datos filtrada (kalman) para generar los pronósticos del modelo, los cuales incluyen:
- Niveles (logaritmos)
- Sumas móviles de de 4T los Productos interno y externo.
- Tasas de variación intertrimestral anualizadas.
- Tasas de variación interanual.

PostProcessing:
Genera los cálculos posteriores que no es posible hacer dentro de la estructura del modelo para un subgrupo de variables:
- Desestacionalización (X12)
- Calculo de brechas y tendencias (HP, lambda=1600)

Posterior a los procesos de tratamiento de datos y simulación se generan las descomposiciones de shocks, gráficas y presentación correspondientes.

Nota: Derivado de la estructura del modelo no fue posible utilizar todas las funciones de SimTools en el proceso de corrimiento.

Departamento de Investigaciones Económicas - Julio 2024. MJGM/JGOR
