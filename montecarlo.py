import numpy as np
import matplotlib.pyplot as plt

# Función para simular el lanzamiento de dos dados
def lanzar_dados():
    return np.random.randint(1, 7) + np.random.randint(1, 7)

# Función mejorada de simulación Monte Carlo
def simulacion_montecarlo_mejorada(n_simulaciones):
    saldo = 0
    suma_anterior = None
    historial_saldos = []
    
    for i in range(n_simulaciones):
        suma_actual = lanzar_dados()

        # Aplicar las reglas de ganancia y pérdida
        if suma_actual in [7, 8, 9, 10, 11, 12]:
            saldo += 10
        elif suma_actual in [2, 3, 4, 5, 6]:
            saldo -= 8
        
        # Aplicar la pérdida por repetición de suma
        if suma_actual == suma_anterior:
            saldo -= 15

        # Actualizar el historial de saldos
        historial_saldos.append(saldo)
        
        # Actualizar la suma anterior
        suma_anterior = suma_actual
    
    return historial_saldos

# Parámetros
n_simulaciones = 30

# Ejecutar la simulación
historial_saldos = simulacion_montecarlo_mejorada(n_simulaciones)

# Mostrar el resultado
print("Historial de saldos acumulados:", historial_saldos)

# Graficar el resultado
plt.plot(historial_saldos)
plt.title("Saldo acumulado en la simulación Monte Carlo")
plt.xlabel("Número de simulaciones")
plt.ylabel("Saldo acumulado")
plt.show()
