import gmpy2 #Para el número primo
from secrets import randbelow  #Para el número aleatorio
import math
import os
import sys
import copy
from contextlib import ExitStack #Librería para abrir varios ficheros a la vez
from glob import glob #Para hacer el ls -l P*.txt

def leer_bloque(f, k, filename):
    def a_int(s):

        v = int.from_bytes(s, byteorder=sys.byteorder)
        v = v.to_bytes((v.bit_length()+7 )//8, byteorder =sys.byteorder )

        fin, inicio = ('', '1') if (len(s) - len(v)) == 0 else ( str(len(s) - len(v)), str(len(str(len(s) - len(v)))+1))

        return int(inicio  +str(int.from_bytes(s, byteorder=sys.byteorder)) + fin)

    def datos_fichero():
        #Leo los datos del fichero, el tamaño del bloque
        tam_tot = os.stat(filename).st_size
        tam_blq = 64
        piezas = tam_tot//tam_blq
        tam_cmp = (tam_blq-4) // (k-1)
        tam_resto = tam_blq - tam_cmp *(k-1)


    piezas, tam_cmp, tam_resto = datos_fichero()
    pos = 0

    while True:
        if pos == (piezas +1):
            break
        pos += 1
        componentes = []

        #Leemos los bloques completos
        for i in range(k):

            s = b''
            if i == k-1:
                s= f.read(tam_resto)
            else:
                s = f.read(tam_cmp)
            componentes.append(a_int(s))

        primo =  gmpy2.next_prime(max(componentes))

        yield componentes, primo


def leer_sistema(dir):

    directorio = dir + "/P?.txt"
    nombre_ficheros = glob(directorio) #Para sacar todos los ficheros con ese nombre

    with ExitStack() as stack:
        ficheros_abiertos = [stack.enter_context(open(fichero)) for fichero in nombre_ficheros] #Para tener un vector con todos los lectores

        while True:

            AB = []
            p = -1
            linea = ""
            v =  []
            for f in ficheros_abiertos:
                linea = f.readline()
                if linea == '':
                    break
                linea = linea[:-1]#quito el '\n'
                linea = linea.split(" ")
                p = int(linea[0])
                v = [int(e) for e in linea[1:]]
                AB.append(v)
            for eq in AB:
                lenSis= len(eq[:-1]) #Cuento las variables
                if lenSis > (len(AB)):
                    print("No se puede encontrar el secreto.")
                    exit(-1)
            yield AB, p
        for j in ficheros_abiertos:
            j.close()


def resolver_sistema(matriz_sistema, primo):#Metodo del montante

    def cofactores(p , i ,j, p_ant):
        cofac = matriz_sistema[p][p]*matriz_sistema[i][j] - matriz_sistema[i][p]* matriz_sistema[p][j]
        cofac = (cofac *p_ant) % primo
        return cofac

    def inverso(B, A): #a el GF(a) B el numero
    #Algoritmode euclides extendido
        g_ant, g_act, u_ant, u_act, v_ant, v_act, i =  A, B, 1,0,0,1,1
        y_f = None

        v = []
        v.append(v_ant)
        v.append(v_act)
        while g_act != 0:
            y_f = g_ant // g_act
            g_f = g_ant - y_f*g_act
            u_f = u_ant - y_f*u_act
            v_f = v_ant -y_f*v_act
            i+=1
            #Actualizamos las variables
            g_ant = g_act
            g_act = g_f

            u_ant = u_act
            u_act = u_f

            v_ant = v_act
            v_act = v_f
            v.append(v_act)
        if v[i-1] < 0:
            v[i-1] += A
        return v[i-1]


    #Ahora resolvemos el sistema
    p_ant = 1
    for i in range(len(matriz_sistema)):
        while matriz_sistema[i][i] == 0:#Comprobación de si el sistema es compatible determinado
            aux = matriz_sistema[i]
            if i+1 >= len(matriz_sistema):
                return None
            matriz_sistema[i]  =matriz_sistema[i +1]
            matriz_sistema[i+1] = aux

        aux_matriz_sistema = [copy.copy(matriz_sistema[i]) for i in range(len(matriz_sistema))]

        for j in range(len(matriz_sistema)):
            if i != j:
                aux_matriz_sistema[j][i] = 0
        for k in range(len(matriz_sistema)):
            for l in range(len(matriz_sistema[i])):
                if k != i :
                    aux_matriz_sistema[k][l]= cofactores(i, k, l, p_ant)
        p_ant = inverso(matriz_sistema[i][i], primo)
        matriz_sistema = aux_matriz_sistema

    pivote =inverso(matriz_sistema[0][0], primo)

    solucion_sistema = []
    div = pivote
    for i in range(len(matriz_sistema)):
        #res = inverso((matriz_sistema[i][-1]*div) %primo , primo)
        solucion_sistema.append((primo+(matriz_sistema[i][-1]*div)) %primo)

    return solucion_sistema


def escribir_secreto(f , secretos):
    for secreto in secretos:
        cadena = str(secreto)
        inicio = int(cadena[0])

        if len(cadena) < 2:
            v= b''

        else:
            fin = 0 if inicio == 1 else  int(cadena[-inicio+1:])

            secreto = int(cadena[1:]) if fin == 0 else int(cadena[1:-inicio+1])

            v= secreto.to_bytes((secreto.bit_length() +7 + 8*fin)//8 , byteorder =sys.byteorder)

        f.write(v)
    return
def construir_hiperplano(S,p):

    def crear_hiperplano(pto_secreto):

        a_i= [ (-p + randbelow(2*p)) for j in range(len(pto_secreto)-1)] #Vector con los coeficientes {a_0, a_2, ..., a_k-2}

        #calculamos el coeficiente independiente, b
        c_i = int(pto_secreto[-1])
        for j in range(len(a_i)):
            c_i -= round(a_i[j]*int(pto_secreto[j]))%p

        #Elhiperplano como un vector
        hiperplano = a_i
        hiperplano.append(-1)  #a_k-1
        hiperplano.append(round((-1)*c_i))

        return hiperplano

    hiperplano = crear_hiperplano(S)

    return hiperplano


def comprobar_datos(k , N, filename):
    #Compruebo si los datos introducidos son buenos
    if k< 2:
        print("No es posible este umbral")
        exit(-1)
    if k > N:
        print("Los datos introducidos no son válidos")
        exit(-1)
    #Abro el fichero con el secreto
    try:
        f = open(filename , 'rb')
    except FileNotFoundError:
        print("El fichero que desea crear los secretos no se ha podido abrir.")
        exit(-2)
    return f


if __name__ == '__main__':
    if len(sys.argv) != 3 and len(sys.argv) != 5:
        print("El número de argumentos no es válido")
        exit(-1)

    if sys.argv[1]== '-c':#Opción de crear secreto

        filename =  sys.argv[2] #Nombre del fichero con el secreto

        #Lo primero voy a decir cual es N y el umbral
        N = int(sys.argv[3])
        k = int(sys.argv[4]) #Determina la dimensión el espacio

        f = comprobar_datos(k, N, filename)

        with ExitStack() as stack:
            #Creo/Vacío los ficheros por si había algo de antes
            archivos =["P" + str(i) + ".txt" for i in range(N)]
            ficheros_abiertos = [stack.enter_context(open(fichero, 'w')) for fichero in archivos]

            #Escribir las ecuaciones en cada fichero
            for (pto_secreto, p) in leer_bloque(f, k, filename):
                for g in ficheros_abiertos:
                    ecuacion = construir_hiperplano( pto_secreto, p)
                    cadena = str(p)+' '
                    for x_i in ecuacion[:-1]:
                        cadena += str(x_i)+' '
                    cadena += str(ecuacion[-1])+"\n"
                    g.write(cadena)

            for escri in ficheros_abiertos:
                escri.close()
            f.close()

    elif sys.argv[1] == '-r': #Opción de reconstruir secreto

        f = open('secreto_reconstruido', "wb")
        dir = sys.argv[2]

        for (matriz_sistema, p) in leer_sistema(dir):

            if len(matriz_sistema) == 0:
                break

            n_min= len(matriz_sistema[0])-1
            #aux_matriz_sistema= copy.copy(matriz_sistema)
            solve= []
            secreto = None
            aux = 0

            while secreto is None:

                aux_matriz_sistema = copy.copy(matriz_sistema)

                for i in range(n_min):
                    j = randbelow(len(aux_matriz_sistema))
                    solve.append(aux_matriz_sistema.pop(j))
                aux += 1
                secreto = resolver_sistema(solve, p)

                if aux > 100:
                    print("No hay secretos suficientes")
                    exit(-1)

            escribir_secreto(f , secreto)
        f.close()
    else:
        print("Argumento no válido,  sólo puede ser -r, para descifrar, y -c, cifrar.")
        exit(-1)
