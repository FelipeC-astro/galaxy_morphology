import matplotlib.pyplot as plt
import numpy as np
from coordinates import converter_ra_dec, angular_distance

def plotar_histograma(dataset, intervalo=250):
    plt.figure(figsize=(6, 4))
    plt.hist(dataset['cz (km/s)'], bins=np.arange(dataset['cz (km/s)'].min(), dataset['cz (km/s)'].max() + intervalo, intervalo))
    plt.xlabel('cz (km/s)')
    plt.ylabel('Frequência')
    plt.title('Histograma de cz (km/s)')
    plt.show()

def plotar_mapa_galaxias(galaxias_e, galaxias_l, centro_ra, centro_dec):
    ra_e, dec_e = converter_ra_dec(galaxias_e['_RAJ2000 (h:m:s)'], galaxias_e['_DEJ2000 (d:m:s)'])
    ra_l, dec_l = converter_ra_dec(galaxias_l['_RAJ2000 (h:m:s)'], galaxias_l['_DEJ2000 (d:m:s)'])
    
    plt.figure(figsize=(8, 6))
    plt.scatter(ra_e, dec_e, marker='*', s=30, color='orange', label='Galáxias Early-Type')
    plt.scatter(ra_l, dec_l, marker='+', s=30, color='purple', label='Galáxias Late-Type')
    plt.scatter(centro_ra, centro_dec, marker='o', s=50, color='black', label='Centro do Aglomerado')
    plt.xlabel('RA (graus)')
    plt.ylabel('DEC (graus)')
    plt.title('Posições das galáxias do aglomerado de Coma')
    plt.grid(True)
    plt.legend(loc='lower left')
    plt.show()

def plotar_histograma_fracoes(galaxias_e, galaxias_l, ra_centro_deg, dec_centro_deg):
    ra_e, dec_e = converter_ra_dec(galaxias_e['_RAJ2000 (h:m:s)'], galaxias_e['_DEJ2000 (d:m:s)'])
    ra_l, dec_l = converter_ra_dec(galaxias_l['_RAJ2000 (h:m:s)'], galaxias_l['_DEJ2000 (d:m:s)'])
    
    distancias_e = angular_distance(ra_e, dec_e, ra_centro_deg, dec_centro_deg)
    distancias_l = angular_distance(ra_l, dec_l, ra_centro_deg, dec_centro_deg)
    
    dist_e = np.radians(distancias_e) * 99
    dist_l = np.radians(distancias_l) * 99
    
    intervalos = np.linspace(min(dist_e.min(), dist_l.min()), max(dist_e.max(), dist_l.max()), num=9)
    
    hist_e, _ = np.histogram(dist_e, bins=intervalos)
    hist_l, _ = np.histogram(dist_l, bins=intervalos)
    
    total_por_bin = hist_e + hist_l
    frac_e = hist_e / total_por_bin
    frac_l = hist_l / total_por_bin
    
    plt.figure(figsize=(8, 6))
    plt.bar(intervalos[:-1], frac_e, width=0.7, color='purple', alpha=0.8, label='Fração de Galáxias Early-Type', align='edge')
    plt.bar(intervalos[:-1], frac_l, width=0.7, color='orange', alpha=0.8, bottom=frac_e, label='Fração de Galáxias Late-Type', align='edge')
    plt.xlabel('Distância (Mpc)')
    plt.ylabel('Fração de galáxias')
    plt.legend(loc='lower left')
    plt.xticks(intervalos.round(2))
    plt.show()
