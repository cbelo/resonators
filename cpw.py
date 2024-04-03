import numpy as np
import scipy.constants as pyc
from scipy.special import ellipk

def fun_k(width, spacing, thickness):
    '''
    Returns the k parameter for the elliptic integrals for the CPW geometry with metallic with, etch spacing and thickness
    '''
    numerator = np.sinh(np.pi*width/(4*thickness))
    denominator = np.sinh(np.pi*(width+2*spacing)/(4*thickness))
    return numerator/denominator

def cpw_cap_diel(k, epsilon_r):
    '''
    Returns the capacitance per unit length of a CPW with dielectric layer with relative permittivity epsilon_r
    '''
    k_prime = np.sqrt(1-k**2)
    K = ellipk(k)
    K_prime = ellipk(k_prime)
    cap = 2*pyc.epsilon_0*(epsilon_r)*K/K_prime
    return cap

def cpw_cap_total(width, spacing, thickness_air, thickness_subs, epsilon_r):
    '''
    Returns the total capacitance per unit length of a CPW with air and dielectric layers
    '''
    k_air = fun_k(width, spacing, thickness_air)#fun_k(width, spacing1, spacing2, thickness_air)
    k_subs = fun_k(width, spacing, thickness_subs)#fun_k(width, spacing1, spacing2, thickness_subs)
    return cpw_cap_diel(k_air, epsilon_r=1) + cpw_cap_diel(k_subs, epsilon_r=epsilon_r)

def cpw_ind_geo(k):
    '''
    Returns the geometric inductance per unit length of a CPW
    '''
    k_prime = np.sqrt(1-k**2)
    K = ellipk(k)
    K_prime = ellipk(k_prime)
    ind = (pyc.mu_0/4)*(K_prime/K)

    return ind

def cpw_ind_geo_total(width, spacing, thickness_air, thickness_subs):
    '''
    Returns the total geometric inductance per unit length of a CPW with air and dielectric layers
    '''
    k_air = fun_k(width, spacing, thickness_air)#fun_k(width_wire, spacing1, thickness_air)
    return cpw_ind_geo(k_air) 

def cpw_ind_kin(ind_kin_sq, width):
    '''
    Returns the kinetic inductance per unit length of a CPW
    '''
    return ind_kin_sq/width

def cpw_ind_total(width, spacing, thickness_air, thickness_subs, ind_kin_sq):
    '''
    Returns the total inductance per unit length of a CPW with air and dielectric layers
    '''
    geometric = cpw_ind_geo_total(width, spacing, thickness_air, thickness_subs)
    kinetic = cpw_ind_kin(ind_kin_sq, width)
    return geometric + kinetic

def impedance_cpw(width, spacing, thickness_air, thickness_subs, epsilon_r, ind_kin_sq):
    '''
    Returns the characteristic impedance of a CPW with air and dielectric layers
    '''
    cap = cpw_cap_total(width, spacing, thickness_air, thickness_subs, epsilon_r)
    ind = cpw_ind_total(width, spacing, thickness_air, thickness_subs, ind_kin_sq)
    return np.sqrt(ind/cap)


def epsilon_eff(width, spacing, thickness_air, thickness_subs, epsilon_r):
    '''
    Returns the effective relative permittivity of a CPW with air and dielectric layers
    '''
    k_air = fun_k(width, spacing, thickness_air)#fun_k(width, spacing1, spacing2, thickness_air)
    K = ellipk(k_air)
    K_prime = ellipk(np.sqrt(1-k_air**2))
    return ((30*pyc.pi /impedance_cpw(width, spacing, thickness_air, thickness_subs, epsilon_r, 0))*(K_prime/K))**2

def f0_cpw(width, spacing, thickness_air, thickness_subs, epsilon_r, ind_kin_sq, length):
    '''
    Returns the resonant frequency of a CPW with air and dielectric layers
    '''
    C = cpw_cap_total(width, spacing, thickness_air, thickness_subs, epsilon_r)
    L = cpw_ind_total(width, spacing, thickness_air, thickness_subs, ind_kin_sq)
    epsilon = epsilon_eff(width, spacing, thickness_air, thickness_subs, epsilon_r)
    f0_LC = 1/(2*np.pi*np.sqrt(L*C)*length)
    f0 = pyc.c/np.sqrt(epsilon)*(1/(2*length))
    print(1/(2*np.sqrt(L*C)*length)*1e-9)
    return np.array([f0, np.pi*f0_LC])

if __name__ == '__main__':
    width = 15.3e-6
    spacing = 3e-6
    thickness_air = 10e-3
    thickness_subs = 500e-6
    epsilon_r = 11.9
    ind_kin_sq = 0#1.2e-6
    length = 1.04e-3
    print('Impedance:', impedance_cpw(width, spacing, thickness_air, thickness_subs, epsilon_r, ind_kin_sq))
    print('Effective relative permittivity:', epsilon_eff(width, spacing, thickness_air, thickness_subs, epsilon_r))
    print('Inductance:', cpw_ind_total(width, spacing, thickness_air, thickness_subs, ind_kin_sq))
    print('Capacitance:', cpw_cap_total(width, spacing, thickness_air, thickness_subs, epsilon_r))
    print('Resonant frequency:', 1e-9*f0_cpw(width, spacing, thickness_air, thickness_subs, epsilon_r, ind_kin_sq, length))