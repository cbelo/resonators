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

def cap_coupling(width_cap, length_horizontal, distance_to_feedline,epsilon_r, thickness_subs):
    '''
    Returns the coupling capacitance of the Schuster rsonator to the feedline
    '''
    k0 = fun_k(width_cap, distance_to_feedline, thickness_subs, simple=True)
    k1 = fun_k(width_cap, distance_to_feedline, thickness_subs, simple=False)
    k_prime0 = np.sqrt(1-k0**2)
    k_prime1 = np.sqrt(1-k1**2)
    K0 = ellipk(k0)
    K0_prime = ellipk(k_prime0)
    K1 = ellipk(k1)
    K1_prime = ellipk(k_prime1)
    cap = 2*pyc.epsilon_0*(epsilon_r-1)*K1/K1_prime +4*pyc.epsilon_0*K0/K0_prime
    return cap*length_horizontal/2

def cap_coupling_new(width_cap, length_horizontal, distance_to_feedline,epsilon_r, thickness_subs, width_feedline):
    
    #Define simpler variables
    s = distance_to_feedline
    w1 = width_cap
    w2 = width_feedline
    h = thickness_subs

    #Calculate the k parameters   
    k0 = np.sqrt(s*(w1+w2+s)/((s+w1)*(s+w2)))
    k0_prime = np.sqrt(1-k0**2)
    K0 = ellipk(k0)
    K0_prime = ellipk(k0_prime)
    num = (np.exp(2*np.pi*(w1+s)/h) -np.exp(2*np.pi*w1/h)) * (np.exp(2*np.pi*(w1+w2+s)/h) -1)
    den = (np.exp(2*np.pi*(w1+w2+s)/h) -np.exp(2*np.pi*w1/h)) * (np.exp(2*np.pi*(w1+s)/h) -1)
    k1 = np.sqrt(num/den)
    k1_prime = np.sqrt(1-k1**2)
    K1 = ellipk(k1)
    K1_prime = ellipk(k1_prime)

    #Calculate the capacitance
    C0 = 2*pyc.epsilon_0*K0_prime/K0
    C1 = pyc.epsilon_0*(epsilon_r-1)*K1_prime/K1
    return (C0 + C1)*length_horizontal

def cap_ground_new(width_cap, length_vertical, distance_to_ground,epsilon_r, thickness_subs):
    #Define simpler variables
    s = distance_to_ground
    w = width_cap
    h = thickness_subs

    #Calculate the k parameters   
    k0 = np.sqrt(w/(w+s))
    k0_prime = np.sqrt(1-k0**2)
    K0 = ellipk(k0)
    K0_prime = ellipk(k0_prime)
    num = np.exp(np.pi*w/h) -1
    den = np.exp(np.pi*(w+s)/h) -1
    k1 = np.sqrt(num/den)
    k1_prime = np.sqrt(1-k1**2)
    K1 = ellipk(k1)
    K1_prime = ellipk(k1_prime)

    #Calculate the capacitance
    C0 = 2*pyc.epsilon_0*K0_prime/K0
    C1 = pyc.epsilon_0*(epsilon_r-1)*K1_prime/K1
    # print('Length',length_vertical)

    return (C0 + C1)*length_vertical*2 #The factor of 2 is because the ground plane is on both sides of the resonator

def cap_ground(width_cap, length_vertical, distance_to_ground,epsilon_r, thickness_subs):
    '''
    Returns the coupling capacitance of the Schuster rsonator to the ground plane
    '''
    k0 = fun_k(width_cap, distance_to_ground, thickness_subs, simple=True)
    k1 = fun_k(width_cap, distance_to_ground, thickness_subs, simple=False)
    k_prime0 = np.sqrt(1-k0**2)
    k_prime1 = np.sqrt(1-k1**2)
    K0 = ellipk(k0)
    K0_prime = ellipk(k_prime0)
    K1 = ellipk(k1)
    K1_prime = ellipk(k_prime1)
    cap = 2*pyc.epsilon_0*(epsilon_r-1)*K1/K1_prime +4*pyc.epsilon_0*K0/K0_prime
    return cap*length_vertical #The factor of 2 is because the ground plane is on both sides of the resonator

def ind_ground_geo(width_ind, length_ind, distance_to_ground):
    '''
    Returns the coupling inductance of the Schuster rsonator to the ground plane
    '''
    k = fun_k(width_ind, distance_to_ground, thickness = 0, simple=True)
    k_prime = np.sqrt(1-k**2)
    K = ellipk(k)
    K_prime = ellipk(k_prime)
    ind = (pyc.mu_0/4)*(K_prime/K)
    #Extra length due to the short to ground
    length_extra = (width_ind + 2*distance_to_ground)/8
    print(length_extra)
    return ind*(length_ind + length_extra) 

def ind_ground_total(width_ind, length_ind, distance_to_ground, ind_kin_sq):
    '''
    Returns the total inductance per unit length of a CPW with air and dielectric layers
    '''
    geometric = ind_ground_geo(width_ind, length_ind, distance_to_ground)
    kinetic = ind_kin_sq/width_ind
    return geometric + kinetic*length_ind

def impedance_Schuster(width_cap, length_horizontal, distance_to_feedline, width_ind, length_ind, distance_to_ground_cap, epsilon_r, thickness_subs, ind_kin_sq):
    '''
    Returns the impedance of the Schuster resonator
    '''
    cap_coupling_val = cap_coupling(width_cap, length_horizontal, distance_to_feedline,epsilon_r, thickness_subs)
    cap_ground_val = cap_ground(width_cap, length_ind, distance_to_ground_cap,epsilon_r, thickness_subs)
    
    distance_to_ground_ind = distance_to_ground_cap + length_horizontal/2+ width_cap/2
    ind_ground_val = ind_ground_total(width_ind, length_ind, distance_to_ground_ind, ind_kin_sq)
    imp0 = np.sqrt(ind_ground_val / (cap_coupling_val + cap_ground_val))
    return imp0

def impedance_Schuster_new(width_cap, length_horizontal_cap, length_vertical_cap, distance_to_feedline, width_ind, length_ind, distance_to_ground_cap, epsilon_r, thickness_subs, ind_kin_sq, width_feedline):
    '''
    Returns the impedance of the Schuster resonator
    '''
    cap_coupling_val = cap_coupling_new(width_cap, length_horizontal_cap, distance_to_feedline,epsilon_r, thickness_subs, width_feedline)
    cap_ground_val = cap_ground_new(width_cap, length_vertical_cap, distance_to_ground_cap,epsilon_r, thickness_subs)
    
    distance_to_ground_ind = distance_to_ground_cap + length_horizontal_cap/2+ width_cap/2
    ind_ground_val = ind_ground_total(width_ind, length_ind, distance_to_ground_ind, ind_kin_sq)
    imp0 = np.sqrt(ind_ground_val / (cap_coupling_val + cap_ground_val))
    return imp0

def resonance_freq_Schuster(width_cap, length_horizontal, distance_to_feedline, width_ind, length_ind, distance_to_ground_cap, epsilon_r, thickness_subs, ind_kin_sq):
    '''
    Returns the resonance frequency of the Schuster resonator
    '''
    cap_coupling_val = cap_coupling(width_cap, length_horizontal, distance_to_feedline,epsilon_r, thickness_subs)
    cap_ground_val = cap_ground(width_cap, length_ind, distance_to_ground_cap,epsilon_r, thickness_subs)
    distance_to_ground_ind = distance_to_ground_cap + length_horizontal/2 + width_cap/2
    print(length_ind)
    ind_ground_val = ind_ground_total(width_ind, length_ind, distance_to_ground_ind, ind_kin_sq)
    print(f'{ind_ground_val*1e9} nH')
    print(f'{1e15*cap_ground_val} fF')
    print(f'{1e15*cap_coupling_val} fF')
    f0 = 1/(2*np.pi*np.sqrt(ind_ground_val*(cap_coupling_val + cap_ground_val)))
    return f0

def resonance_freq_Schuster_new(width_cap, length_horizontal_cap, length_vertical_cap, distance_to_feedline, width_ind, length_ind, distance_to_ground_cap, epsilon_r, thickness_subs, ind_kin_sq, width_feedline):
    '''
    Returns the resonance frequency of the Schuster resonator
    '''
    cap_coupling_val = cap_coupling_new(width_cap, length_horizontal_cap, distance_to_feedline,epsilon_r, thickness_subs, width_feedline)
    cap_ground_val = cap_ground_new(width_cap, length_vertical_cap, distance_to_ground_cap,epsilon_r, thickness_subs)
    distance_to_ground_ind = distance_to_ground_cap + length_horizontal_cap/2 + width_cap/2
    ind_ground_val = ind_ground_total(width_ind, length_ind, distance_to_ground_ind, ind_kin_sq)
    print(f'{ind_ground_val*1e9} nH')
    print(f'{1e15*cap_ground_val} fF')
    print(f'{1e15*cap_coupling_val} fF')
    f0 = 1/(2*np.pi*np.sqrt(ind_ground_val*(cap_coupling_val + cap_ground_val)))
    return f0

if __name__ == '__main__':
    width_cap = 25e-6
    length_horizontal = 500e-6
    distance_to_feedline = 3e-6
    width_ind = 1e-6
    length_ind = 200e-6
    distance_to_ground_cap = 3e-6
    epsilon_r = 11.9
    thickness_subs = 500e-6
    ind_kin_sq = 100e-12
    print(f' Imp: {impedance_Schuster(width_cap, length_horizontal, distance_to_feedline, width_ind, length_ind, distance_to_ground_cap, epsilon_r, thickness_subs, ind_kin_sq)} Ohm')
    print(f'Freq: {1e-9*resonance_freq_Schuster(width_cap, length_horizontal, distance_to_feedline, width_ind, length_ind, distance_to_ground_cap, epsilon_r, thickness_subs, ind_kin_sq)} GHz')