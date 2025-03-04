import pandas as pd

def carregar_dados(caminho):
    colunas = ["_RAJ2000 (h:m:s)", "_DEJ2000 (d:m:s)", "Num", "CGCG", "NGC_IC", "UGC", "cz (km/s)", 
               "RA1950 (h:m:s)", "DE1950 (d:m:s)", "B25.5(mag)", "SB25.5(mag)", "5 (arcsec)", "pe", 
               "_RA.icrs (h:m:s)", "_DE.icrs (d:m:s)"]
    
    dataset = pd.read_csv(caminho, header=None, comment='#', names=colunas, sep='|', usecols=range(1, len(colunas) + 1))
    dataset['cz (km/s)'] = dataset['cz (km/s)'].str.strip().replace('', '0').astype(float)
    return dataset[dataset['cz (km/s)'] != 0]

def filtrar_galaxias(dataset, v_aglomerado=6925, dispersao=1000):
    limite_inf, limite_sup = v_aglomerado - 2 * dispersao, v_aglomerado + 2 * dispersao
    return dataset[(dataset['cz (km/s)'] >= limite_inf) & (dataset['cz (km/s)'] <= limite_sup)]

