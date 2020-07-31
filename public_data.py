steel_bars = {'Tag': ['D[mm]', 'A[mm²]'],
              '#3': [9.5, 71],
              '#4': [12.7, 129],
              '#5': [15.9, 199],
              '#6': [19.1, 284],
              '#7': [22.2, 387],
              '#8': [25.4, 510],
              '#9': [28.7, 645],
              '#10': [32.3, 819]}

rebar_methods = {'Intercalda': 'intercalado',
                 'Esquinero': 'esquinero'}

column_design_factors = {'defConc': 0.003,  # Ɛ  deformación unitaria de falla del concreto
                         'deftrac': 0.005,  # Ɛ  deformación unitaria para falla por tracción
                         'ficomp': 0.65,    # ɸc para sección controlada por compresión
                         'fitracc': 0.90,   # ɸt para sección controlada por tracción
                         'alfa': 1.0}    # α  factor de amplificación por sobreresistencia del acero

iteration_limit = {'ia': 100,  # iteraciones por ángulo
                   'ic': 50,   # iteraciones por profundidad c
                   'ip': 20,   # iteraciones por ángulo y c
                   'errC': 0.0001,  # error máximo para profundidad c
                   'errA': 0.01,    # error máximo para el ángulo
                   'errP': 0.0005   # error máximo para fuerza P = errP*Ag*f'c
                   }
