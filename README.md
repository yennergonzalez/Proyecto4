
### Modelos Probabilísticos de Señales y Sistemas - IE0405


#### Yenner Josué González Araya - B83375
# Proyecto 4 - Modulación Digital IQ
Este proyecto consistió en tres partes: la simulación de la transmisión de una imagen bajo un esquema de
modulación 16-QAM, la realización de pruebas estacionaridad y ergodicidad sobre la señal modulada y
la determinación y graficación de la densidad espectral de potencia de la señal modulada 'senal_Tx'.

## Parte 1 - Modulación 16-QAM
Para esto se modificó el código ya implementado en el enunciado del proyecto para expandir
la funcionalidad del código de una modulación BPSK a 16-QAM. 
Esto se logró añadiendo una segunda señal portadora en la función 'modulador_QAM' de manera
que se tuviera una señal 'I' correspondiente a coseno y una señal 'Q' asociada a una forma de
onda seno. Luego se fueron recorriendo los bits obtenidos al realizar la conversión RGB a bits
de la imagen, y de acuerdo a los pares de bits se define la amplitud que tendrá cada señal.

Luego, se suman las señales 'I' y 'Q' para agruparlas en una única señal transmitida 'señal_Tx'.
Esta función moduladora retorna esta señal transmitida al igual que la potencia de la señal y la 
señal portadora.

Para la demodulación se utilizó el mismo método utilizado para la modulación BPSK, por lo que
cualquier 'pseudo-energía' mayor a cero en un período de la portadora se traduce como un '1'
y los otros casos se traducen como un '0'. Finalmente esta función retorna la secuencia de
bits recibidos y la señal demodulada.

A continuación se presentan las funciones de modulación y demodulación implementadas para la simulación
de la transmisión con modulación 16-QAM:
  
### Modulación:
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

### Demodulación:

    import numpy as np

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

Para la simulación realizada se obtuvieron los siguientes resultados:
Formas de onda de la señal modulada, la señal modulada con ruido y la señal demodulada, respectivamente:
![ScreenShot](https://raw.githubusercontent.com/yennergonzalez/Proyecto4/main/Formas_de_Onda_16QAM.png)

Comparación de la imagen transmitida y recuperada 
![ScreenShot](https://raw.githubusercontent.com/yennergonzalez/Proyecto4/main/Imagen_Enviada_y_Recuperada_16QAM.png)

## Parte 2 - Estacionaridad y ergodicidad
 
Para esta parte se debía analizar la señal modulada 'senal_Tx' y  realizar un análisis de ergodicidad y estacionaridad. 

La estacionaridad ocurre cuando el promedio estadístico es igual al promedio temporal, por lo que se
calcularon ambos para la señal 'senal_Tx' mediante el siguiente código:
    
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

Para este código se obtuvo la siguiente salida:

    El promedio temporal de la  señal Tx es: 0.004128658
    El promedio estadístico de la señal Tx es: 0.004128835.
    El porcentaje de error entre los promedios estadístico y temporal de la señal Tx es del: 0.000017685 %.
    Es claro que ambos promedios son sumamente similares entre sí, por lo que se puede concluir que la señal Tx es ergódica, y dado que la ergodicidad es también una forma de estacionaridad, la señal Tx es estacionaria también.

Donde se puede observar que los valores obtenidos para el promedio temporal y el promedio estadístico son prácticamente iguales, por lo que se puede concluir que la señal 'senal_Tx' es ergódica y dado que la ergodicidad es una forma de estacionaridad se puede concluir también que esta señal es  estacionaria también.

## Parte 3 - Densidad espectral de potencia

Para la tercera parte se debía obtener una gráfica de la densidad espectral de potencia 'senal_Tx'.
Esta se obtuvo mediante la ejecución del siguiente código:

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
Para este se obtuvo como resultado la siguiente gráfica:
![ScreenShot](https://raw.githubusercontent.com/yennergonzalez/Proyecto4/main/Densidad_Espectral_de_Potencia_16QAM.png)
