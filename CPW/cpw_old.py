import numpy as np
import scipy.constants as pyc

#### Lukas functions
def fun_k_Lukas(width, spacing):
    return width / (width + 2 * spacing)

def cpw_cap_diel_Lukas(k, epsilon_r):
    if k <= 1 / np.sqrt(2):
        kprime = np.sqrt(1 - k ** 2)
        f = np.pi / np.log(2 * (1 + np.sqrt(kprime)) / (1 - np.sqrt(kprime)))
    else:
        f = np.log(2 * (1 + np.sqrt(k)) / (1 - np.sqrt(k))) / np.pi
    return 2 * pyc.epsilon_0 * epsilon_r * f

def cpw_cap_total_Lukas(k, epsilon_r):
    return cpw_cap_diel_Lukas(k, epsilon_r=1) + cpw_cap_diel_Lukas(k, epsilon_r=epsilon_r)

def cpw_ind_geo_Lukas(k):
    if k <= 1 / np.sqrt(2):
        kprime = np.sqrt(1 - k ** 2)
        f = np.pi / np.log(2 * (1 + np.sqrt(kprime)) / (1 - np.sqrt(kprime)))
    else:
        f = np.log(2 * (1 + np.sqrt(k)) / (1 - np.sqrt(k))) / np.pi
    return pyc.mu_0 / 4 / f

def cpw_ind_kin_Lukas(ind_kin_sq, width): ## zero order approximation                          
    return ind_kin_sq / width

def cpw_ind_total_Lukas(k, ind_kin_sq, width):
    return cpw_ind_geo_Lukas(k) + cpw_ind_kin_Lukas(ind_kin_sq, width)

def cpw_bare_resonator_Lukas(width, spacing, length, ind_kin_sq, epsilon_r):
    k = fun_k(width, spacing)
    cap_m = cpw_cap_total_Lukas(k, epsilon_r)
    ind_m = cpw_ind_total_Lukas(k, ind_kin_sq, width)
    freq = 1/(4*length*np.sqrt(ind_m*cap_m))
    imp = np.sqrt(ind_m/cap_m)
    return np.array([freq, imp], dtype=object)

def cpw_resonator_length_Lukas(width, spacing, freq_nw, ind_kin_sq, ind_nw, epsilon_r):
    # properties of unloaded cpw resonator                          
    k = fun_k(width, spacing)
    cap_m = cpw_cap_total_Lukas(k, epsilon_r)
    print('Capacitance per meter: {} F/m'.format(cap_m))
    ind_m = cpw_ind_total_Lukas(k, ind_kin_sq, width)
    imp0 = np.sqrt(ind_m / cap_m)
    print(f'imp0 = {imp0} Ohm')
    # length = 1 / (2*np.pi * freq_nw * np.sqrt(ind_m * cap_m))
    # based on characteristic impedance and target frequency, find the right resonator length                          
    f0 = freq_nw / (1 - 4 * freq_nw * ind_nw / imp0)
    length = 1 / (4 * f0 * np.sqrt(ind_m * cap_m))
    return length

#### Miguel functions (based on Simons' book)
'''
Metallic layer with thickness 'thickness', central conductor with width 'width' and spacing 'spacing1' adn 'spacing2'
In the usual CPW configuration, spacing1 = spacing2
Above the metal, there's air with thickness 'thickness_air' 
Below the metal, a dielectric layer with thickness 'thickness_subs' and relative permittivity 'epsilon_r
'''

from scipy.special import ellipk

def fun_k(width, spacing1, spacing2, thickness):
    numerator = np.tanh(np.pi*width/(4*thickness))
    denominator = np.tanh(np.pi*(width+spacing1+spacing2)/(4*thickness))
    return numerator/denominator

def fun_k_new(width, spacing, thickness):
    return np.sinh(np.pi*width/(4*thickness))/np.sinh(np.pi*(width+2*spacing)/(4*thickness))

def cpw_cap_diel(k, epsilon_r):
    k_prime = np.sqrt(1-k**2)
    K = ellipk(k)
    K_prime = ellipk(k_prime)
    cap = 2*pyc.epsilon_0*(epsilon_r-1)*K/K_prime
    return cap

def cpw_cap_total(width, spacing1, spacing2, thickness_air, thickness_subs, epsilon_r):
    k_air = fun_k_new(width, spacing1, thickness_air)#fun_k(width, spacing1, spacing2, thickness_air)
    k_subs = fun_k_new(width, spacing1, thickness_subs)#fun_k(width, spacing1, spacing2, thickness_subs)
    return cpw_cap_diel(k_air, epsilon_r=1) + cpw_cap_diel(k_subs, epsilon_r=epsilon_r)

def cpw_ind_geo(k):
    k_prime = np.sqrt(1-k**2)
    K = ellipk(k)
    K_prime = ellipk(k_prime)
    ind = (pyc.mu_0/4)*(K_prime/K)
    return ind

def cpw_ind_geo_total(width_wire, spacing1, spacing2, thickness_air, thickness_subs):
    k_air = fun_k_new(width_wire, spacing1, thickness_air)#fun_k(width_wire, spacing1, spacing2, thickness_air)
    k_subs = fun_k_new(width_wire, spacing1, thickness_subs)#fun_k(width_wire, spacing1, spacing2, thickness_subs)
    return cpw_ind_geo(k_air) + cpw_ind_geo(k_subs)

def cpw_ind_kin(ind_kin_sq, width_wire):
    return ind_kin_sq/width_wire

def cpw_ind_total(width_wire, spacing1, spacing2, thickness_air, thickness_subs, ind_kin_sq):
    geometric = cpw_ind_geo_total(width_wire, spacing1, spacing2, thickness_air, thickness_subs)
    kinetic = cpw_ind_kin(ind_kin_sq, width_wire)
    return geometric + kinetic

if __name__ == '__main__':
    width = 15.3e-6
    spacing1 = 3e-6
    spacing2 = spacing1
    thickness_air = 1e-3
    thickness_subs = 500e-6
    epsilon_r = 11.9
    ind_kin_sq = 0#4e-12
    length = 1e-3
    L = cpw_ind_total(width, spacing1, spacing2, thickness_air, thickness_subs, ind_kin_sq)
    print(f'Inductance = {L*1e9} nH/m')
    C = cpw_cap_total(width, spacing1, spacing2, thickness_air, thickness_subs, epsilon_r)
    Z = np.sqrt(L/C)
    print(f'Characteristic impedance = {Z} Ohm')
    print(f'Capacitance = {C*1e15} fF/m')
    k = fun_k_Lukas(width, spacing1)
    C_L = cpw_cap_total_Lukas(k, epsilon_r)
    L_L = cpw_ind_total_Lukas(k, ind_kin_sq, width)
    print(f'Capacitance = {C_L*1e15} fF/m')
    print(f'Inductance = {L_L*1e9} nH/m')
    Z_L = np.sqrt(L_L/C_L)
    print(f'Characteristic impedance = {Z_L} Ohm')
    # freq = 1/(2*np.pi*np.sqrt(L*C))
    # print(f'Frequency = {freq*1e-9} GHz')