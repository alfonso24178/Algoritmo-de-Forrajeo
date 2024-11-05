from bacteria import bacteria
from chemiotaxis import chemiotaxis
import numpy
import random  # Importar la librería random

# Inicialización de variables
poblacion = []
path = r"C:\Users\herna\OneDrive\Escritorio\Escritorio\Willars\multifasta.fasta"
numeroDeBacterias = 5
numRandomBacteria = 1
iteraciones = 30
tumbo = 1                                              # Número de gaps a insertar
nado = 3
chemio = chemiotaxis()
veryBest = bacteria(path)                # Mejor bacteria   
tempBacteria = bacteria(path)            # Bacteria temporal para validaciones
original = bacteria(path)                # Bacteria original sin gaps
globalNFE = 0      # Número de evaluaciones de la función objetivo
dAttr = 0.1 # 0.1
wAttr = 0.2 # 0.2
hRep = dAttr
wRep = 10    # 10
no_improvement_cycles = 0  # Inicializa el contador de ciclos sin mejora
previous_best_fitness = veryBest.fitness  # Guardar el fitness inicial

# Función para dispersar bacterias
def disperse_bacteria(bacteria_population, num_eliminate):
    # Eliminar bacterias aleatorias de la población
    return random.sample(bacteria_population, len(bacteria_population) - num_eliminate)

# Función para clonar la mejor bacteria
def clonaBest(veryBest, best):
    veryBest.matrix.seqs = numpy.array(best.matrix.seqs)
    veryBest.blosumScore = best.blosumScore
    veryBest.fitness = best.fitness
    veryBest.interaction = best.interaction

# Función para validar las secuencias
def validaSecuencias(path, veryBest):
    # Clonar a veryBest en tempBacteria   
    tempBacteria.matrix.seqs = numpy.array(veryBest.matrix.seqs)
    # Descartar los gaps de cada secuencia
    for i in range(len(tempBacteria.matrix.seqs)):
        tempBacteria.matrix.seqs[i] = tempBacteria.matrix.seqs[i].replace("-", "")
    # Validar que las secuencias originales sean iguales a las secuencias de tempBacteria
    for i in range(len(tempBacteria.matrix.seqs)):
        if tempBacteria.matrix.seqs[i] != original.matrix.seqs[i]:
            print("*****************Secuencias no coinciden********************")
            return

# Parámetros adaptativos
def update_attract_repel(iteration, max_iter):
    attract = max(1 - (iteration / max_iter), 0.1)  # decrece con las iteraciones
    repel = min((iteration / max_iter), 0.9)  # aumenta con las iteraciones
    return attract, repel

# Función de eliminación adaptativa
def adaptive_elimination(bacteria_population, fitness, iteration, no_improvement_cycles, threshold=0.01):
    if no_improvement_cycles >= threshold:
        # Eliminar más bacterias si no hay mejora significativa
        num_eliminate = int(0.5 * len(bacteria_population))  # eliminar 50% de las bacterias
        bacteria_population = disperse_bacteria(bacteria_population, num_eliminate)
    return bacteria_population

# Población inicial de bacterias
for i in range(numeroDeBacterias):
    poblacion.append(bacteria(path))

# Bucle principal de iteraciones
for iteration in range(iteraciones):
    # Actualizar los valores de atracción y repulsión
    dAttr, hRep = update_attract_repel(iteration, iteraciones)
    
    # Procesar cada bacteria en la población
    for bacteria in poblacion:
        bacteria.tumboNado(tumbo)
        bacteria.autoEvalua()

    # Realizar el proceso de chemio-taxis
    chemio.doChemioTaxis(poblacion, dAttr, wAttr, hRep, wRep)
    globalNFE += chemio.parcialNFE

    # Encontrar la mejor bacteria
    best = max(poblacion, key=lambda x: x.fitness)
    
    # Verificar si ha habido mejora en el fitness
    if best.fitness > previous_best_fitness:
        no_improvement_cycles = 0  # Reiniciar el contador si hay mejora
        clonaBest(veryBest, best)
        previous_best_fitness = best.fitness
    else:
        no_improvement_cycles += 1  # Incrementar el contador si no hay mejora

    # Aplicar la eliminación adaptativa si no hay mejoras
    poblacion = adaptive_elimination(poblacion, best.fitness, iteration, no_improvement_cycles)
    
    # Mostrar resultados
    print("Interacción: ", veryBest.interaction, "Fitness: ", veryBest.fitness, " NFE:", globalNFE)

    # Eliminar bacterias y clonar nuevas
    chemio.eliminarClonar(path, poblacion)
    chemio.insertRamdomBacterias(path, numRandomBacteria, poblacion)  # Insertar bacterias aleatorias
    print("Población: ", len(poblacion))

# Mostrar el genoma de la mejor bacteria
veryBest.showGenome()
validaSecuencias(path, veryBest)



