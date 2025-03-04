import pandas as pd
from data_processing import carregar_dados, filtrar_galaxias
from coordinates import converter_ra_dec
from plotting import plotar_histograma, plotar_mapa_galaxias, plotar_histograma_fracoes

# Carregar os dados
dataset = carregar_dados("./MODIFICADOCatalogComaCluster.txt")
plotar_histograma(dataset)

# Filtrar galáxias
galaxias_filtradas = filtrar_galaxias(dataset)
galaxias_e = galaxias_filtradas[galaxias_filtradas['pe'] == 'E']
galaxias_l = galaxias_filtradas[galaxias_filtradas['pe'] == 'L']

# Centro do aglomerado
ra_centro, dec_centro = converter_ra_dec(pd.Series(['12 59 48.7']), pd.Series(['+27 58 50']))

# Plotar os gráficos
plotar_mapa_galaxias(galaxias_e, galaxias_l, ra_centro[0], dec_centro[0])
plotar_histograma_fracoes(galaxias_e, galaxias_l, ra_centro[0], dec_centro[0])
