# Proyecto 4 - Modelos Probabilísticos de Señales y Sistemas
# Yenner Josué González Araya - B83375

# Funciones del enunciado del proyecto:

from PIL import Image
import numpy as np

def fuente_info(imagen):
    '''Una función que simula una fuente de
    información al importar una imagen y 
    retornar un vector de NumPy con las 
    dimensiones de la imagen, incluidos los
    canales RGB: alto x largo x 3 canales

    :param imagen: Una imagen en formato JPG
    :return: un vector de pixeles
    '''
    img = Image.open(imagen)
    
    return np.array(img)

def rgb_a_bit(array_imagen):
    '''Convierte los pixeles de base 
    decimal (de 0 a 255) a binaria 
    (de 00000000 a 11111111).

    :param imagen: array de una imagen 
    :return: Un vector de (1 x k) bits 'int'
    '''
    # Obtener las dimensiones de la imagen
    x, y, z = array_imagen.shape
    
    # Número total de elementos (pixeles x canales)
    n_elementos = x * y * z

    # Convertir la imagen a un vector unidimensional de n_elementos
    pixeles = np.reshape(array_imagen, n_elementos)

    # Convertir los canales a base 2
    bits = [format(pixel, '08b') for pixel in pixeles]
    bits_Rx = np.array(list(''.join(bits)))
    
    return bits_Rx.astype(int)

def canal_ruidoso(senal_Tx, Pm, SNR):
    '''Un bloque que simula un medio de trans-
    misión no ideal (ruidoso) empleando ruido
    AWGN. Pide por parámetro un vector con la
    señal provieniente de un modulador y un
    valor en decibelios para la relación señal
    a ruido.

    :param senal_Tx: El vector del modulador
    :param Pm: Potencia de la señal modulada
    :param SNR: Relación señal-a-ruido en dB
    :return: La señal modulada al dejar el canal
    '''
    # Potencia del ruido generado por el canal
    Pn = Pm / pow(10, SNR/10)

    # Generando ruido auditivo blanco gaussiano (potencia = varianza)
    ruido = np.random.normal(0, np.sqrt(Pn), senal_Tx.shape)

    # Señal distorsionada por el canal ruidoso
    senal_Rx = senal_Tx + ruido

    return senal_Rx

def bits_a_rgb(bits_Rx, dimensiones):
    '''Un blque que decodifica el los bits
    recuperados en el proceso de demodulación

    :param: Un vector de bits 1 x k 
    :param dimensiones: Tupla con dimensiones de la img.
    :return: Un array con los pixeles reconstruidos
    '''
    # Cantidad de bits
    N = len(bits_Rx)

    # Se reconstruyen los canales RGB
    bits = np.split(bits_Rx, N / 8)

    # Se decofican los canales:
    canales = [int(''.join(map(str, canal)), 2) for canal in bits]
    pixeles = np.reshape(canales, dimensiones)

    return pixeles.astype(np.uint8)

# --------------------------------------------------------------------

# Modulación 16 - QAM

# Definición de la función de  modulación en 16-QAM

def modulador_QAM(bits, fc, mpp):
    # Método para simular un esquema de modulación digital QAM-16.

    # 1. Parámetros señal de información (bits)
    N = len(bits) # Cantidad de bits  

    # 2. Construcción señales portadoras c(t)
    Tc = 1/fc   # Período [s]
    t_periodo = np.linspace(0, Tc, mpp) # Muestras por período

    # Inicialización de las señales portadoras 
    portadora_I = np.cos(2*np.pi*fc*t_periodo)
    portadora_Q = np.sin(2*np.pi*fc*t_periodo)

    # 3. Inicialización señales moduladas s(t)
    t_simulacion = np.linspace(0, N*Tc, N*mpp) 
    senal_Tx_I = np.zeros(t_simulacion.shape)
    senal_Tx_Q = np.zeros(t_simulacion.shape)
    
    #4. Prueba inicial de indexación de bits
    for i in range(0, len(bits), 2): # Amplitudes
        
        # A1 --- portadora I
        if bits[i]==0 and bits[i+1]==0:
            senal_Tx_I[i*mpp : (i+1)*mpp]=-3*portadora_I

        if bits[i]==0 and bits[i+1]==1:
            senal_Tx_I[i*mpp : (i+1)*mpp]=-1*portadora_I

        if bits[i]==1 and bits[i+1]==0:
            senal_Tx_I[i*mpp : (i+1)*mpp]=3*portadora_I

        if bits[i]==1 and bits[i+1]==1:
            senal_Tx_I[i*mpp : (i+1)*mpp]=1*portadora_I

        # A2 --- portadora Q
        if bits[i]==0 and bits[i+1]==0:
            senal_Tx_Q[(i+1)*mpp : (i+2)*mpp]=3*portadora_Q

        if bits[i]==0 and bits[i+1]==1:
            senal_Tx_Q[(i+1)*mpp : (i+2)*mpp]=1*portadora_Q

        if bits[i]==1 and bits[i+1]==0:            
            senal_Tx_Q[(i+1)*mpp : (i+2)*mpp]=-3*portadora_Q

        if bits[i]==1 and bits[i+1]==1:
            senal_Tx_Q[(i+1)*mpp : (i+2)*mpp]=1*portadora_Q

     # Potencias de cada señal modulada
    Pm_I = (1 / (N*Tc)) * np.trapz(pow(senal_Tx_I, 2), t_simulacion)
    Pm_Q = (1 / (N*Tc)) * np.trapz(pow(senal_Tx_Q, 2), t_simulacion)

    Pm = Pm_I+Pm_Q
    senal_Tx = senal_Tx_Q + senal_Tx_I
    portadora = portadora_Q + portadora_I    
    
    return senal_Tx, Pm, portadora

def demodulador_QAM(senal_Rx, portadora, mpp):
    '''Un método que simula un bloque demodulador
    de señales, bajo un esquema 16-QAM. El criterio
    de demodulación se basa en decodificación por 
    detección de energía.

    :param senal_Rx: La señal recibida del canal
    :param portadora: La onda portadora c(t)
    :param mpp: Número de muestras por periodo
    :return: Los bits de la señal demodulada
    '''
    # Cantidad de muestras en senal_Rx
    M = len(senal_Rx)

    # Cantidad de bits (símbolos) en transmisión
    N = int(M / mpp)

    # Vector para bits obtenidos por la demodulación
    bits_Rx = np.zeros(N)

    # Vector para la señal demodulada
    senal_demodulada = np.zeros(senal_Rx.shape)

    # Pseudo-energía de un período de la portadora
    Es = np.sum(portadora**2)
    
    # Demodulación
    for i in range(N):
        # Producto interno de dos funciones
        producto = senal_Rx[i*mpp : (i+1)*mpp] * portadora
        Ep = np.sum(producto) 
        senal_demodulada[i*mpp : (i+1)*mpp] = producto

        # Criterio de decisión por detección de energía
        if Ep > 0:
            bits_Rx[i] = 1
        else:
            bits_Rx[i] = 0
    
    return bits_Rx.astype(int), senal_demodulada

# Transmisión de una imagen bajo modulación 16-QAM

import numpy as np
import matplotlib.pyplot as plt
import time

# Parámetros
fc = 5000  # frecuencia de la portadora
mpp = 20   # muestras por periodo de la portadora
SNR = 5   # relación señal-a-ruido del canal

# Iniciar medición del tiempo de simulación
inicio = time.time()

# 1. Importar y convertir la imagen a trasmitir
imagen_Tx = fuente_info('arenal.jpg')
dimensiones = imagen_Tx.shape

# 2. Codificar los pixeles de la imagen
bits_Tx = rgb_a_bit(imagen_Tx)

# 3. Modular la cadena de bits usando el esquema 16-QAM
senal_Tx, Pm, portadora = modulador_QAM(bits_Tx, fc, mpp)

# 4. Se simula el paso por el canal ruidoso
senal_Rx = canal_ruidoso(senal_Tx, Pm, SNR)

# 5. Demodular los bits usando el esquema 16-QAM
bits_Rx, senal_demodulada = demodulador_QAM(senal_Rx, portadora, mpp)

# 6. Pasar de bits a RGB
imagen_Tx_QAM = bits_a_rgb(bits_Tx, dimensiones)
Fig = plt.figure(figsize=(10,6))

imagen_Rx_QAM = bits_a_rgb(bits_Rx, dimensiones)
Fig = plt.figure(figsize=(10,6))

# Mostrar imagen transmitida
ax = Fig.add_subplot(1, 2, 1)
imgplot = plt.imshow(imagen_Tx_QAM)
ax.set_title('Transmitido')

# Mostrar imagen recuperada
ax = Fig.add_subplot(1, 2, 2)
imgplot = plt.imshow(imagen_Rx_QAM)
ax.set_title('Recuperado')
Fig.tight_layout()

plt.imshow(imagen_Rx_QAM)

# 7. Formas de onda

# Visualizar el cambio entre las señales
fig, (ax2, ax3, ax4) = plt.subplots(nrows=3, sharex=True, figsize=(14, 7))

# La señal modulada por 16-QAM
ax2.plot(senal_Tx[0:600], color='g', lw=2) 
ax2.set_ylabel('$s(t)$')

# La señal modulada al dejar el canal
ax3.plot(senal_Rx[0:600], color='b', lw=2) 
ax3.set_ylabel('$s(t) + n(t)$')

# La señal demodulada
ax4.plot(senal_demodulada[0:600], color='m', lw=2) 
ax4.set_ylabel('$b^{\prime}(t)$')
ax4.set_xlabel('$t$ / milisegundos')
fig.tight_layout()
plt.show()


# 8. Calcular número de errores
errores = sum(abs(bits_Tx - bits_Rx))
BER = errores/len(bits_Tx)
print('{} errores, para un BER de {:0.4f}.'.format(errores, BER))


# Parte 2 - Estacionaridad y ergodicidad

# Verificación de estacionaridad y ergodicidad

# Para verificar ergodicidad se deben calcular los promedios temporales y los estadísticos, y si son iguales o muy parecidos, el proceso se denomina ergódico.

# Cálculo promedio temporal.

# Muestras de la señal
Muestras = len(senal_Tx)
# Número de símbolos
Simbolos = Muestras // mpp
# Tiempo del periodo
T_Periodo = 1 / fc
# Periodo en segundos de un ciclo
t_simulacion = np.linspace(0, Simbolos*T_Periodo, Simbolos*mpp)
T_integral = t_simulacion[len(t_simulacion)-1]
# Se integra de 0 a infinito por lo que se divide por T y no por T/2.
med_temp = (1 / T_integral) * np.trapz(senal_Tx, t_simulacion)
print("El promedio temporal de la  señal Tx es: {:0.9f}".format(med_temp))


# Cálculo promedio estadístico.

med_estad = np.mean(senal_Tx)
print('El promedio estadístico de la señal Tx es: {:0.9f}.'.format(med_estad))

# Resultados y conclusión.
porcentaje_error =  abs(med_estad-med_temp)*100
print('El porcentaje de error entre los promedios estadístico y temporal de la señal Tx es del: {:0.9f} %.'.format(porcentaje_error))
print('Es claro que ambos promedios son sumamente similares entre sí, por lo que se puede concluir que la señal Tx es ergódica, y dado que la ergodicidad es también una forma de estacionaridad, la señal Tx es estacionaria también.')



# Parte 3 - Densidad espectral de potencia

from scipy import fft

# Transformada de Fourier
senal_f = fft(senal_Tx)

# Muestras de la señal
Nm = len(senal_Tx)

# Número de símbolos (198 x 89 x 8 x 3)
Ns = Nm // mpp

# Tiempo del símbolo = periodo de la onda portadora
Tc = 1 / fc

# Tiempo entre muestras (período de muestreo)
Tm = Tc / mpp

# Tiempo de la simulación
T = Ns * Tc

# Espacio de frecuencias
f = np.linspace(0.0, 1.0/(2.0*Tm), Nm//2)

# Gráfica
plt.plot(f, 2.0/Nm * np.power(np.abs(senal_f[0:Nm//2]), 2))
plt.xlim(0, 10000)
plt.grid()
plt.title("Densidad espectral de potencia de la señal Tx")
plt.xlabel("Frecuencia (Hz)")
plt.ylabel("Amplitud")
plt.show()