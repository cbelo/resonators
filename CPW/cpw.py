import numpy as np
import scipy.constants as pyc
from scipy.special import ellipk

def fun_k(width, spacing, thickness, simple = False):
    '''
    Returns the k parameter for the elliptic integrals for the CPW geometry with metallic with, etch spacing and thickness
    '''
    if simple:
        return width/(width+2*spacing)
    else:
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

def cpw_cap_total(width, spacing, thickness_air, thickness_subs, epsilon_r,  simple =False):
    '''
    Returns the total capacitance per unit length of a CPW with air and dielectric layers
    '''
    k_air = fun_k(width, spacing, thickness_air, simple = simple)#fun_k(width, spacing1, spacing2, thickness_air)
    k_subs = fun_k(width, spacing, thickness_subs, simple = simple)#fun_k(width, spacing1, spacing2, thickness_subs)
    return cpw_cap_diel(k_air, epsilon_r=2) + cpw_cap_diel(k_subs, epsilon_r=epsilon_r-1)

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
    vph = 1/np.sqrt(C*L)
    f0 = vph/(2*length)
    return f0

def cpw_simpler(width, spacing, thickness_air, thickness_subs, epsilon_r, ind_kin_sq):
    '''
    Returns the characteristic impedance of a CPW with air and dielectric layers
    '''
    k0 = fun_k(width, spacing, thickness_subs, simple = True)
    k1 = fun_k(width, spacing, thickness_subs, simple = True)
    k_prime0 = np.sqrt(1-k0**2)
    k_prime1 = np.sqrt(1-k1**2)
    K0 = ellipk(k0)
    K0_prime = ellipk(k_prime0)
    K1 = ellipk(k1)
    K1_prime = ellipk(k_prime1)
    C1 = 2*pyc.epsilon_0*(epsilon_r-1)*K1/K1_prime 
    C_air = 4*pyc.epsilon_0*(K0/K0_prime)
    C_CPW = C1 + C_air
    epsilon = C_CPW/C_air
    Z = 1/(pyc.c*C_air*np.sqrt(epsilon))
    L = (pyc.mu_0/4)*(K0_prime/K0)
    
    return C_CPW, epsilon, Z, L
                   


if __name__ == '__main__':
    width = 100e-6
    spacing = 50e-6
    thickness_air = 10e-3
    thickness_subs = 500e-6
    epsilon_r = 11.9
    ind_kin_sq = 0e-12
    length = 1.04e-3
    print('Impedance:', impedance_cpw(width, spacing, thickness_air, thickness_subs, epsilon_r, ind_kin_sq))
    print('Effective relative permittivity:', epsilon_eff(width, spacing, thickness_air, thickness_subs, epsilon_r))
    print('Inductance (nH/m):', 1e9*cpw_ind_total(width, spacing, thickness_air, thickness_subs, ind_kin_sq))
    print('Capacitance (pF/m):', 1e12*cpw_cap_total(width, spacing, thickness_air, thickness_subs, epsilon_r, simple =False))
    print('Resonant frequency (GHz):', 1e-9*f0_cpw(width, spacing, thickness_air, thickness_subs, epsilon_r, ind_kin_sq, length))
    print('\n Simplified model:')
    C_CPW, epsilon, Z, L = cpw_simpler(width, spacing, thickness_air, thickness_subs, epsilon_r, ind_kin_sq)
    print('Impedance:', Z)
    print('Effective relative permittivity:', epsilon)
    print('Inductance (nH/m):', 1e9*L)
    print('Capacitance (pF/m):', 1e12*C_CPW)