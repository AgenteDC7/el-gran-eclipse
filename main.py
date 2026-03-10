
# *****************************************************************
#    PRACTICA 1 - ECLIPSES - PARADIGMAS DE PROGRAMACIÓN 2025-26
#    PLANTILLA DE LA PRÁCTICA
#    Apellidos, Nombre: Clermontel Villalba, Daniel
#    Apellidos, Nombre: Gómez Blanco, Marvin Ayelén
# *****************************************************************

# ephem: https://rhodesmill.org/pyephem/index.html
# ephem: https://pypi.org/project/ephem/
# ephem: pip install ephem 

# global_land_mask: https://pypi.org/project/global-land-mask/
# global_land_mask: pip install global-land-mask
# global_land_mask: globe.is_land(lat, lon), lat-lon in degrees

# ansi: https://gist.github.com/fnky/458719343aabd01cfb17a3a4f7296797

import math
from ephem import *
from global_land_mask import globe
from time import sleep, time

# ************************ PARÁMETROS DEL PROGRAMA *****************************

N_SEG = 7                      # Número de segmentos en los que se divide cada intervalo para buscar el mínimo
PRECISION = (1.0/86400)/N_SEG  # Precisión (en fracciones de dia) en el cálculo de la fecha del eclipse (1 segundo)
SEP_LIM = 1.5 * degree         # Separación mínima para que pueda existir un eclipse de sol (radianes)
LIM_LAT = 70                   # Latitud máxima que se representa (grados)
N_FIL = 50                     # Número de filas de las matrices de mapa y animación
N_COL = 128                    # Número de columnas de las matrices de mapa y animación
ANIM_TAM_Y = 0.03              # Tamaño vertical (altura en radianes) del cuadro de animación
ANIM_TAM_X = ANIM_TAM_Y*N_COL/N_FIL  # Tamaño horizontal (azimut en radianes) del cuadro de animación

# COLORES ANSI-256 PARA MAPA (función col_mapa)

COLS_MAPA = ((18, 25, 32, 75),    # Gradiente agua
             (235, 88, 130, 178))  # Gradiente tierra
COL_SITIO = 196                   # Marca del sitio (rojo intenso)

# COLORES ANSI-256 PARA ANIMACIÓN (función col_anim)

COLS_ANIM = (75, 226, 242, 242, 130, 214, 172, 172)


# ************************ FUNCIONES PREDEFINIDAS *****************************

def cls():
    """ Borrado de pantalla (códigos ANSI) """
    print("\x1b[H\x1b[2J\x1b[3J")


def traduce_latlon(txt: str) -> (str | float, str | float):
    """ Traduce una posición geográfica en formato Google Maps (41°28'40.2"N 4°35'53.8"W) a una
        tupla de strings (latitud, longitud) con el formato de ephem
    :param txt: La posición obtenida de Google Maps
    :return: Una tupla de strings o floats (latitud, longitud)
    """
    if ',' in txt:
        return tuple(float(s) * degree for s in txt.split(", "))
    lat, lon = txt.split(' ')
    lat = lat.replace('°', ':').replace("'", ':')
    lon = lon.replace('°', ':').replace("'", ':')
    return (lat[0:-2] if lat[-1] == 'N' else '-' + lat[0:-2],
            lon[0:-2] if lon[-1] == 'E' else '-' + lon[0:-2])


def dist_ang(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    """  Calcula la distancia angular entre dos puntos en una esfera
    :param lon1: Longitud o Ascensión Recta o Azimut del primer punto (radianes)
    :param lat1: Latitud o Declinación o Altura del primer punto (radianes)
    :param lon2: Longitud o Ascensión Recta o Azimut del segundo punto (radianes)
    :param lat2: Latitud o Declinación o Altura del segundo punto (radianes)
    :return: Distancia angular (radianes) entre los puntos
    """
    return acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))


def ocultacion(r1: float, r2: float, d: float) -> float:
    """ Calcula el porcentaje de area ocultada de un círculo por otro círculo
    :param r1: Radio del círculo ocultado
    :param r2: Radio del círculo ocultador
    :param d: Distancia entre los centros de ambos círculos
    :return: Un valor real entre 0 (no hay solapamiento) y 100 (totalmente ocultado)
    """
    if d > r1 + r2:  # Sin solapamiento
        return 0.0
    if d + r2 <= r1:  # Ocultador mas pequeño y completamente dentro del ocultado
        return 100 * r2 * r2 / (r1 * r1)
    if d + r1 <= r2:  # Ocultador mayor y tapa completamente al ocultado
        return 100.0
    # Area de la region de interseccion
    ai = r1 * r1 * acos((d * d + r1 * r1 - r2 * r2) / (2 * d * r1)) + \
        r2 * r2 * acos((d * d - r1 * r1 + r2 * r2) / (2 * d * r2)) - \
        0.5 * math.sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
    # Ratio del area del primer círculo
    return 100 * ai / (pi * r1 * r1)


def sep_alt_ocult(fec: Date, obj1: Body, obj2: Body, obs: Observer | None = None) -> (float, float, float):
    """ Calcula la distancia angular entre los cuerpos en esa fecha, la altura sobre el horizonte del primer cuerpo y
    el porcentaje de ocultación del primer cuerpo por parte del segundo.
    Si no se proporciona punto de observación los resultados son geocéntricos.
    :param fec: Fecha
    :param obj1: Cuerpo celeste ocultado
    :param obj2: Cuerpo celeste ocultador
    :param obs: Lugar de observación (None si el cálculo es geocéntrico)
    :return: Una tupla con la separación angular entre los cuerpos (radianes), altura sobre el horizonte (radianes) y
     porcentaje de ocultación (0-100)
    """
    if obs is None:
        # Cálculo geocéntrico, usamos RA/DEC y devolvemos una altura positiva cualquiera
        obj1.compute(fec)
        obj2.compute(fec)
        sep, alt = dist_ang(obj1.g_ra, obj1.g_dec, obj2.g_ra, obj2.g_dec), 0.00001
    else:
        # Cálculo topocéntrico, usamos Azimut/Altura
        obs.date = fec
        obj1.compute(obs)
        obj2.compute(obs)
        sep, alt = dist_ang(obj1.az, obj1.alt, obj2.az, obj2.alt), obj1.alt
    return sep, alt, ocultacion(obj1.radius, obj2.radius, sep)


def en_circulo(cx: float, cy: float, r: float, x: float, y: float):
    """ Comprueba si el punto (x,y) se encuentra en el círculo de centro (cx,cy) y radio r
    :param cx: Centro (x) del círculo
    :param cy: Centro (y) del círculo
    :param r: Radio del círculo
    :param x: Posición (x) del punto
    :param y: Posición (y) del punto
    :return: True si el punto está dentro del círculo
    """
    return (cx - x) * (cx - x) + (cy - y) * (cy - y) <= r * r


def fec_local(fec: Date) -> str:
    """ Convierte una fecha a texto en hora local """
    f = fec if isinstance(fec, Date) else Date(fec)
    return localtime(f).isoformat(' ', 'seconds') if f >= 25567.5 else str(f)


def interv(cen: float, tam: float, n: int) -> [float]:
    """ Devuelve una lista de los puntos centrales de un intervalo divido en segmentos contiguos iguales
    :param cen: Centro del intervalo
    :param tam: Longitud del intervalo
    :param n: Número de segmentos en que se divide el intervalo
    :return: Lista de posiciones de los puntos centrales de cada uno de los n segmentos, ordenada
    """
    x0 = cen - tam * (n - 1) / (2 * n)
    return [x0 + i * tam / n for i in range(n)]


def col_mapa(ocul: float, lat: float, lon: float) -> int:
    """ Devuelve el color de un pixel del mapa dado el nivel de ocultación y de si es tierra u océano
    :param ocul: El porcentaje (0-100) de ocultación
    :param lat: Latitud del punto (radianes)
    :param lon: Longitud del punto (radianes)
    :return: Color ANSI-256 del pixel
    """
    return COLS_MAPA[int(globe.is_land(lat/degree, lon/degree))][int(1.8 * math.log10(102 - ocul))]


def col_anim(en_obj1: bool, en_obj2: bool, bajo_horiz: bool = False) -> int:
    """ Devuelve el color del un pixel de animación según las condiciones del punto
    :param en_obj1: El punto está dentro del círculo que representa el primer cuerpo celeste
    :param en_obj2: El punto está dentro del círculo que representa el segundo cuerpo celeste
    :param bajo_horiz: El punto está bajo el horizonte
    :return: Color ANSI-256 del pixel
    """
    return COLS_ANIM[en_obj1 | en_obj2 << 1 | bajo_horiz << 2]


# ************************ FUNCIONES AUXILIARES DEL ALUMNO *******************************


# ************************ CABECERAS DE FUNCIONES AUXILIARES *****************************

def lista_eclipses(tipo: int, fec_ini: Date, fec_fin: Date, obj1: Body, obj2: Body, obs: Observer, min_ocu: float)\
        -> [(Date, Date, float, float, float, float, str)]:
    """ Devuelve la lista de los eclipses que se pueden dar entre dos fechas, filtrando por separación
    mínima desde el punto de vista geocéntrico y ocultación mínima desde el punto de vista del observador
    :param tipo: Tipo de eclipse (1 -> Eclipse de Sol)
    :param fec_ini: Fecha inicial
    :param fec_fin: Fecha final
    :param obj1: Cuerpo celeste ocultado
    :param obj2: Cuerpo celeste ocultador
    :param obs: Punto de observación
    :param min_ocu: Mínimo valor de ocultación para ser incluido (porcentaje, 1-100)
    :return: Enumeración de datos de eclipses: (fecha geo, fecha top, sep geo, sep top, altitud, ocultación, mensaje),
     las separaciones y altitud en grados.
    """
    pass


def menu_principal(tipo: int, lis: [(Date, Date, float, float, float, float, str)],
                   obj1: Body, obj2: Body, obs: Observer):
    """ Muestra el menu de acciones (listado, mapa y animación) sobre los eclipses
    :param tipo: Tipo de cálculo (0 -> Eclipses de Sol)
    :param lis: Lista de posibles eclipses y su información asociada
    :param obj1: Cuerpo celeste ocultado
    :param obj2: Cuerpo celeste ocultador
    :param obs: Lugar de observación
    """
    pass


# FUNCIÓN PRINCIPAL

def main():
    print("ECLIPSES - PARADIGMAS 2025-26")
    print("1. Eclipse de Sol")
    tipo = input("Escoja tipo de cálculo (1): ")
    tipo = int(tipo) if tipo else 1

    # Datos necesarios
    fec_ini = input("Fecha inicial (formato año-mes-dia): ")
    fec_fin = input("Fecha final (formato año-mes-dia):   ")
    loc_geo = input("Localidad, posicion geografica: ")
    loc_alt = input("Localidad, altura (m): ")
    min_ocu = input("Porcentaje minimo de ocultacion (0..100): ")
    # Traducción y sustitución de valores por defecto
    fec_ini = Date(fec_ini) if fec_ini else now()
    fec_fin = Date(fec_fin) if fec_fin else Date(fec_ini+3653)
    # La localidad por defecto es Valladolid
    lat, lon = traduce_latlon(loc_geo) if loc_geo else (41.66308134*degree, -4.70494676*degree)
    loc_alt = int(loc_alt) if loc_alt else 0
    min_ocu = int(min_ocu) if min_ocu else 0
    # Creación del punto de observación y de los cuerpos celestes
    sitio = Observer()
    sitio.lat = lat
    sitio.lon = lon
    sitio.elevation = loc_alt
    obj1, obj2 = Sun(), Moon()
    # Obtención de la lista de eclipses
    lis = lista_eclipses(tipo, fec_ini, fec_fin, obj1, obj2, sitio, min_ocu)
    # Menu principal del programa
    menu_principal(tipo, lis, obj1, obj2, sitio)


if __name__ == '__main__':
    main()