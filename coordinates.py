import numpy as np
import pandas as pd

def converter_ra_dec(ra, dec):
    def ra_to_deg(ra):
        h, m, s = map(float, ra.split())
        return (h + m / 60 + s / 3600) * 15
    
    def dec_to_deg(dec):
        d, m, s = map(float, dec.split())
        sinal = -1 if d < 0 else 1
        return sinal * (abs(d) + m / 60 + s / 3600)
    
    return ra.apply(ra_to_deg).values, dec.apply(dec_to_deg).values

def angular_distance(ra1, dec1, ra2, dec2):
    ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])
    delta_ra = ra2 - ra1
    delta_dec = dec2 - dec1
    a = np.sin(delta_dec / 2)**2 + np.cos(dec1) * np.cos(dec2) * np.sin(delta_ra / 2)**2
    return 2 * np.arcsin(np.sqrt(a)) * (180 / np.pi)

