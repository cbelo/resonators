import numpy as np
import scipy.constants as pyc
from scipy.special import ellipk, ellipkm1


def fun_k_CPW(width, spacing, thickness, simple = False):
    '''
    Returns the k parameter for the elliptic integrals for the CPW geometry with metallic with, etch spacing and thickness
    If simple is True, it returns the k parameter for the CPW geometry with infinite substrate (thickness = infinity).
    '''
    if simple:
        return width/(width+2*spacing)
    else:######### FIX FOR LOW THIKNESS
        if np.pi*(width+2*spacing)/(4*thickness)>709:
            k = np.exp(-np.pi*(spacing)/(2*thickness))
        else:
            numerator = np.sinh(np.pi*width/(4*thickness))
            denominator = np.sinh(np.pi*(width+2*spacing)/(4*thickness))
            k = numerator/denominator
        return k

def C_air_CPW(width, spacing):
    '''
    Returns the capacitance through air of the CPW to the ground plane.
    '''
    k = fun_k_CPW(width, spacing, thickness = 0, simple=True)
    k_prime = np.sqrt(1-k**2)
    K = ellipk(k)
    K_prime = ellipk(k_prime)
    C_air = 4*pyc.epsilon_0*K/K_prime
    return C_air

def cap_CPW_per_length(width, spacing, epsilon_r, thickness_subs):
    '''
    Returns the  capacitance of the CPW to the ground plane per length.
    Eplison_r is the relative permittivity of the substrate. It can be a number or an array.
    The thickness of the substrate h must have the same dimension than epsilon_r.
    Order of epsilon_r & h: The first element is the layer closer to the metallic layer.
    '''

    if isinstance(epsilon_r, (int, float)):
        epsilon_r = np.array([epsilon_r])
        thickness_subs = np.array([thickness_subs])

    C_air = C_air_CPW(width, spacing)
    cap_contributions = [C_air]

    for i in range(len(epsilon_r)):
        k = fun_k_CPW(width, spacing, thickness_subs[i])
        k_prime = np.sqrt(1-k**2)
        K = ellipk(k)
        K_prime = ellipkm1(k_prime)
        if i == len(epsilon_r)-1:
            e_r = epsilon_r[i]-1
        else:
            e_r = epsilon_r[i]-epsilon_r[i+1]
        cap = 2*pyc.epsilon_0*e_r*K/K_prime
        cap_contributions.append(cap)
    return sum(cap_contributions)

def ind_CPW_per_length(width, spacing, ind_kin_sq):
    '''
    Returns the inductance of the CPW to the ground plane per length.
    Eplison_r is the relative permittivity of the substrate. It can be a number or an array.
    The thickness of the substrate h must have the same dimension than epsilon_r.
    Order of epsilon_r & h: The first element is the layer closer to the metallic layer.
    '''
    k = fun_k_CPW(width, spacing, thickness=0, simple=True)
    k_prime = np.sqrt(1-k**2)
    K = ellipk(k)
    K_prime = ellipk(k_prime)
    ind_geo = (pyc.mu_0/4)*(K_prime/K)
    ind_kin = ind_kin_sq/width
    return ind_kin + ind_geo

def impedance_CPW(width, spacing, epsilon_r, thickness_subs, ind_kin_sq):
    '''
    Returns the impedance of the CPW.
    Eplison_r is the relative permittivity of the substrate. It can be a number or an array.
    The thickness of the substrate h must have the same dimension than epsilon_r.
    Order of epsilon_r & h: The first element is the layer closer to the metallic layer.
    '''
    C = cap_CPW_per_length(width, spacing, epsilon_r, thickness_subs)
    L = ind_CPW_per_length(width, spacing, ind_kin_sq)
    return np.sqrt(L/C)

def resonance_freq_CPW(width, spacing, epsilon_r, thickness_subs, ind_kin_sq, length_CPW):
    '''
    Returns the resonance frequency of the CPW.
    Eplison_r is the relative permittivity of the substrate. It can be a number or an array.
    The thickness of the substrate h must have the same dimension than epsilon_r.
    Order of epsilon_r & h: The first element is the layer closer to the metallic layer.
    '''
    C = cap_CPW_per_length(width, spacing, epsilon_r, thickness_subs)
    L = ind_CPW_per_length(width, spacing, ind_kin_sq)
    vph = 1/np.sqrt(C*L)
    print(length_CPW)
    f0 = vph/(2*np.pi*length_CPW)
    return f0

def epsilon_eff_CPW(width, spacing, epsilon_r, thickness_subs):
    '''
    Returns the effective relative permittivity of the CPW.
    Eplison_r is the relative permittivity of the substrate. It can be a number or an array.
    The thickness of the substrate h must have the same dimension than epsilon_r.
    Order of epsilon_r & h: The first element is the layer closer to the metallic layer.
    '''
    C = cap_CPW_per_length(width, spacing, epsilon_r, thickness_subs)
    C_air = C_air_CPW(width, spacing)
    return C/C_air

if __name__ == '__main__':
    width = 100e-6
    spacing = 50e-6
    thickness_air = 10e-3
    thickness_subs = 500e-6
    epsilon_r = 11.9
    ind_kin_sq = 0e-12
    length = 1.04e-3
    print('Impedance:', impedance_CPW(width, spacing, epsilon_r, thickness_subs, ind_kin_sq))
    print('Freq:', 1e-9*resonance_freq_CPW(width, spacing, epsilon_r, thickness_subs, ind_kin_sq, length))
    print('Epsilon_eff:', epsilon_eff_CPW(width, spacing, epsilon_r, thickness_subs))
    print('Cap per length:', 1e12*cap_CPW_per_length(width, spacing, epsilon_r, thickness_subs))
    print('Ind per length:', 1e9*ind_CPW_per_length(width, spacing, ind_kin_sq))