import numpy as np
import scipy.constants as pyc
from scipy.special import ellipk

def fun_k_CPW(width, spacing, thickness, simple = False):
    '''
    Returns the k parameter for the elliptic integrals for the CPW geometry with metallic with, etch spacing and thickness
    If simple is True, it returns the k parameter for the CPW geometry with infinite substrate (thickness = infinity).
    '''
    if simple:
        return width/(width+2*spacing)
    else:
        numerator = np.sinh(np.pi*width/(4*thickness))
        denominator = np.sinh(np.pi*(width+2*spacing)/(4*thickness))
        return numerator/denominator

def fun_k_asymmetric_stripline(w1, w2, s, h):
    '''
    Computes the argument of the elliptic integrals for the an asymmetric coplanar stripline
    This is used to compute Cc in the Schuster resonator.
    
    '''
    num = (np.exp(2*np.pi*(w1+s)/h) -np.exp(2*np.pi*w1/h)) * (np.exp(2*np.pi*(w1+w2+s)/h) -1)
    den = (np.exp(2*np.pi*(w1+w2+s)/h) -np.exp(2*np.pi*w1/h)) * (np.exp(2*np.pi*(w1+s)/h) -1)
    k = np.sqrt(num/den)
    return k
    
def fun_k_stripline_to_infinite_ground(w,s,h):
    '''
    Computes the argument of the elliptic integrals for the coplanar stripline to infinite ground
    This is used to compute C0 in the Schuster resonator.
    '''
    num = np.exp(np.pi*w/h) -1
    den = np.exp(np.pi*(w+s)/h) -1
    k = np.sqrt(num/den)
    return k

def C_air_asymmetric_stripline(w1, w2, s):
    '''
    Returns the capacitance through air of the Schuster resonator to the feedline.
    '''
    k0 = np.sqrt(s*(w1+w2+s)/((s+w1)*(s+w2)))
    k0_prime = np.sqrt(1-k0**2)
    K0 = ellipk(k0)
    K0_prime = ellipk(k0_prime)
    C_air = 2*pyc.epsilon_0*K0_prime/K0
    return C_air

def C_air_stripline_to_infinite_ground(w, s):
    '''
    Returns the capacitance through air of the Schuster resonator to the ground.
    '''
    k0 = np.sqrt(w/(w+s))
    k0_prime = np.sqrt(1-k0**2)
    K0 = ellipk(k0)
    K0_prime = ellipk(k0_prime)
    C_air = 2*pyc.epsilon_0*K0_prime/K0
    return C_air

def cap_coupling(width_cap, horizontal_length_cap, distance_to_feedline, width_feedline, epsilon_r, thickness_subs):
    '''
    Returns the coupling capacitance of the Schuster resonator to the feedline.
    Eplison_r is the relative permittivity of the substrate. It can be a number or an array.
    The thickness of the substrate h must have the same dimension than epsilon_r.
    Order of epsilon_r & h: The first element is the layer closer to the metallic layer.
    '''

    #Define simpler variables
    s = distance_to_feedline
    w1 = width_cap
    w2 = width_feedline

    if isinstance(epsilon_r, (int, float)):
        epsilon_r = np.array([epsilon_r])
        thickness_subs = np.array([thickness_subs])

    capacitance_contributions = []

    C_air =  C_air_asymmetric_stripline(w1, w2, s)
    capacitance_contributions.append(C_air*horizontal_length_cap)

    for i in range(len(epsilon_r)):
        k = fun_k_asymmetric_stripline(w1, w2, s, thickness_subs[i])
        k_prime = np.sqrt(1-k**2)
        K = ellipk(k)
        K_prime = ellipk(k_prime)
        if i == len(epsilon_r)-1:
            e_r = epsilon_r[i]-1
        else:
            e_r = epsilon_r[i]-epsilon_r[i+1]
        cap = pyc.epsilon_0*e_r*K_prime/K
        capacitance_contributions.append(cap*horizontal_length_cap)

    return sum(capacitance_contributions)

def cap_ground(width_cap, vertical_length_cap, distance_to_ground, epsilon_r, thickness_subs):
    '''
    Returns the coupling capacitance of the Schuster resonator to the ground plane.
    Eplison_r is the relative permittivity of the substrate. It can be a number or an array.
    The thickness of the substrate h must have the same dimension than epsilon_r.
    Order of epsilon_r & h: The first element is the layer closer to the metallic layer.
    '''
    #Define simpler variables
    s = distance_to_ground
    w = width_cap

    if isinstance(epsilon_r, (int, float)):
        epsilon_r = np.array([epsilon_r])
        thickness_subs = np.array([thickness_subs])

    capacitance_contributions = []

    C_air =  C_air_stripline_to_infinite_ground(w, s)
    capacitance_contributions.append(C_air*vertical_length_cap*2) # Factor of 2 since there are two ground planes

    for i in range(len(epsilon_r)):
        k = fun_k_stripline_to_infinite_ground(w, s, thickness_subs[i])
        k_prime = np.sqrt(1-k**2)
        K = ellipk(k)
        K_prime = ellipk(k_prime)
        if i == len(epsilon_r)-1:
            e_r = epsilon_r[i]-1
        else:
            e_r = epsilon_r[i]-epsilon_r[i+1]
        cap = pyc.epsilon_0*e_r*K_prime/K
        capacitance_contributions.append(cap*vertical_length_cap*2)

    return sum(capacitance_contributions)

def ind_ground_geo(width_ind, length_ind, distance_to_ground):
    '''
    Returns the coupling inductance of the Schuster rsonator to the ground plane
    '''
    k = fun_k_CPW(width_ind, distance_to_ground, thickness = 0, simple=True) 
    #test simple and no simple!!!!!!!!
    k_prime = np.sqrt(1-k**2)
    K = ellipk(k)
    K_prime = ellipk(k_prime)
    ind = (pyc.mu_0/4)*(K_prime/K)
    #Extra length due to the short to ground
    length_extra = (width_ind + 2*distance_to_ground)/8
    return ind*(length_ind + length_extra) 

def ind_ground_total(width_ind, length_ind, distance_to_ground, ind_kin_sq):
    '''
    Returns the total inductance per unit length of a CPW with air and dielectric layers
    '''
    geometric = ind_ground_geo(width_ind, length_ind, distance_to_ground)
    kinetic = ind_kin_sq/width_ind
    return geometric + kinetic*length_ind ## Test with length extra!

def impedance_Schuster(width_cap, horizontal_length_cap, distance_to_feedline, width_feedline, #Cc
                         vertical_length_cap, distance_to_ground_cap, #Cg
                         width_ind, length_ind, ind_kin_sq, #Lg
                         epsilon_r, thickness_subs): #Epsilon_r, h
    '''
    Returns the impedance of the Schuster resonator.
    Eplison_r is the relative permittivity of the substrate. It can be a number or an array.
    The thickness of the substrate h must have the same dimension than epsilon_r.
    Order of epsilon_r & h: The first element is the layer closer to the metallic layer.
    '''
    Cc = cap_coupling(width_cap, horizontal_length_cap, distance_to_feedline, width_feedline, epsilon_r, thickness_subs)
    Cg = cap_ground(width_cap, vertical_length_cap, distance_to_ground_cap, epsilon_r, thickness_subs)
    distance_to_ground_ind = distance_to_ground_cap + horizontal_length_cap/2+ width_cap/2
    Lg = ind_ground_total(width_ind, length_ind, distance_to_ground_ind, ind_kin_sq)
    return np.sqrt(Lg/(Cg+Cc))

def resonance_freq_Schuster(width_cap, horizontal_length_cap, distance_to_feedline, width_feedline, #Cc
                         vertical_length_cap, distance_to_ground_cap, #Cg
                         width_ind, length_ind, ind_kin_sq, #Lg
                         epsilon_r, thickness_subs): #Epsilon_r, h
    '''
    Returns the resonance frequency of the Schuster resonator.
    Eplison_r is the relative permittivity of the substrate. It can be a number or an array.
    The thickness of the substrate h must have the same dimension than epsilon_r.
    Order of epsilon_r & h: The first element is the layer closer to the metallic layer.
    '''
    Cc = cap_coupling(width_cap, horizontal_length_cap, distance_to_feedline, width_feedline, epsilon_r, thickness_subs)
    Cg = cap_ground(width_cap, vertical_length_cap, distance_to_ground_cap, epsilon_r, thickness_subs)
    distance_to_ground_ind = distance_to_ground_cap + horizontal_length_cap/2 + width_cap/2
    Lg = ind_ground_total(width_ind, length_ind, distance_to_ground_ind, ind_kin_sq)
    print(f'{Lg*1e9} nH')
    print(f'{1e15*Cg} fF')
    print(f'{1e15*Cc} fF')
    return 1/(2*np.pi*np.sqrt(Lg*(Cc+Cg)))

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
        k = fun_k_CPW(width, spacing, thickness_subs)
        k_prime = np.sqrt(1-k**2)
        K = ellipk(k)
        K_prime = ellipk(k_prime)
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
    f0 = vph/(2*length_CPW)
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
    SeparationTlineResonator = 0
    NumberOfResonators = 1
    FeedlineWidth = 120
    FeedlineLength  = 1000
    FeedlineGap = 3
    FeedlineTaperLength =  50
    BondpadWidth =  300
    BondpadLength = 200
    BondpadGap = 50

    #Resonator parameters
    CapacitorHorizontalLength = np.ones(NumberOfResonators)*120
    CapacitorVerticalLength = np.ones(NumberOfResonators)*200
    CapacitorWidth = np.ones(NumberOfResonators)*15

    NumberOfBends = np.ones(NumberOfResonators, dtype=int)*50
    InductorVerticalLength = np.ones(NumberOfResonators)*11
    InductorHorizontalLength = np.ones(NumberOfResonators)*50
    InductorEndLength = 50*np.ones(NumberOfResonators)
    InductorTotalLength = InductorVerticalLength*(NumberOfBends+4) + InductorEndLength +InductorHorizontalLength*NumberOfBends
    InductorWidth = np.ones(NumberOfResonators)*1
    TaperWidth = np.ones(NumberOfResonators)*20
    SpacingC0 = np.ones(NumberOfResonators)*10
    SpacingCc = np.ones(NumberOfResonators)*(2)
    TaperLength = np.ones(NumberOfResonators)*5

    epsilon_r = 11.9
    thickness_subs = 500e-6
    ind_kin_sq = 100e-12
    # print(f' Imp: {impedance_Schuster(CapacitorWidth*1e-6, CapacitorHorizontalLength*1e-6, (SeparationTlineResonator + FeedlineGap + SpacingCc)*1e-6 , FeedlineWidth*1e-6,
    #                                    CapacitorVerticalLength*1e-6, SpacingC0*1e-6,
    #                                    InductorWidth*1e-6, InductorTotalLength*1e-6, ind_kin_sq,
    #                                    epsilon_r, thickness_subs)} Ohm')
    print(f'Freq: {1e-9*resonance_freq_Schuster(CapacitorWidth*1e-6, CapacitorHorizontalLength*1e-6, (SeparationTlineResonator + FeedlineGap + SpacingCc)*1e-6 , FeedlineWidth*1e-6,
                                       CapacitorVerticalLength*1e-6, SpacingC0*1e-6,
                                       InductorWidth*1e-6, InductorTotalLength*1e-6, ind_kin_sq,
                                       epsilon_r, thickness_subs)} GHz')


    ############# CPW TESTS ################
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
    