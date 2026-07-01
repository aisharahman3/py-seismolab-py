import math

def energy(magnitude):
    return 10 ** (1.5 * magnitude + 4.8)

def moment_magnitude(m0):
    return (2.0 / 3.0) * (math.log10(m0) - 9.1)

def energy_ratio(m1, m2):
    return 10 ** (1.5 * (m1 - m2))

def b_value(magnitudes, mc=None):
    """Aki (1965) maximum-likelihood b-value for a catalog above completeness Mc."""
    mags = [m for m in magnitudes if mc is None or m >= mc]
    if not mags:
        raise ValueError("no events at or above completeness magnitude mc")
    if mc is None:
        mc = min(mags)
    mean_m = sum(mags) / len(mags)
    if mean_m == mc:
        raise ValueError("degenerate catalog: all magnitudes equal mc, b-value undefined")
    return math.log10(math.e) / (mean_m - mc)

def epicentral_distance(lat1, lon1, lat2, lon2, radius=6371.0):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * radius * math.asin(math.sqrt(a))
