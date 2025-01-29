import numpy as np
import scipy.constants as pyc
from scipy.special import ellipk, ellipkm1

'''
This script contains the functions to compute the impedance and resonance frequency of the Schuster resonator, using conformal mapping and elliptic integrals.
The functions are based on the book  'Coplanar Waveguide Circuits, Components, and Systems' by Rainee N. Simons.
These functions are valid for different substrate layers with different relative permittivities and thicknesses.

The functions are:
- fun_k_CPW: Returns the k parameter for the elliptic integrals for the CPW geometry with metallic width, etch spacing and substrate thickness (Section 2.2.3).
- fun_k_asymmetric_stripline: Computes the argument of the elliptic integrals for the an asymmetric coplanar stripline (Section 6.2.3 in the book).
- fun_k_stripline_to_infinite_ground: Computes the argument of the elliptic integrals for the coplanar stripline to infinite ground (Section 6.2.4 in the book).
- C_air_asymmetric_stripline: Returns the capacitance through air of the Schuster resonator to the feedline (Section 6.2.3 in the book).
- C_air_stripline_to_infinite_ground: Returns the capacitance through air of the Schuster resonator to the ground (Section 6.2.4 in the book).
- cap_coupling: Returns the coupling capacitance of the Schuster resonator to the feedline (Section 6.2.3 in the book).
- cap_ground: Returns the coupling capacitance of the Schuster resonator to the ground plane (Section 6.2.4 in the book).
- ind_ground_geo: Returns the coupling inductance of the Schuster resonator to the ground plane (Section 12.4.3 and 9.3 in the book).
- ind_ground_total: Returns the total (geometric+kinetic) inductance per unit length of a CPW with air and dielectric layers.
- impedance_Schuster: Returns the impedance of the Schuster resonator.
- resonance_freq_Schuster: Returns the resonance frequency of the Schuster resonator.

'''

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
    if 2*np.pi*(w1+w2+s)/h > 709:
        k = 1
    else:
        num = np.float64((np.exp(2*np.pi*(w1+s)/h) -np.exp(2*np.pi*w1/h)) * (np.exp(2*np.pi*(w1+w2+s)/h) -1))
        den = np.float64(np.exp(2*np.pi*(w1+w2+s)/h) -np.exp(2*np.pi*w1/h)) * (np.exp(2*np.pi*(w1+s)/h) -1)
        k = np.sqrt(num/den)
    return k
    
def fun_k_stripline_to_infinite_ground(w,s,h):
    '''
    Computes the argument of the elliptic integrals for the coplanar stripline to infinite ground
    This is used to compute C0 in the Schuster resonator.
    '''
    if np.pi*w/h > 709:
        k = np.sqrt(np.exp(-np.pi*s/h))
    else:
        num = np.array(np.exp(np.pi*w/h) -1)
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
    C_air = 2*pyc.epsilon_0*(K0_prime/K0)#**(-1)
    return C_air

def C_air_stripline_to_infinite_ground(w, s):
    '''
    Returns the capacitance through air of the Schuster resonator to the ground.
    '''
    # print('Air distance',s)
    # print('Air width',w)
    ### From simmons book 
    k0 = np.sqrt(w/(w+s))
    k0_prime = np.sqrt(1-k0**2)

    K0 = ellipk(k0)
    K0_prime = ellipk(k0_prime)
    C_air = 2*pyc.epsilon_0*(K0_prime/K0)#**(-1) 
    # print(f' C_air: {C_air*1e12} pF/m')
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
        if np.allclose(k, 1):
            K = ellipkm1(k)
        else:
            K = ellipk(k)
        K_prime = ellipk(k_prime)
        if i == len(epsilon_r)-1:
            e_r = epsilon_r[i]-1
        else:
            e_r = epsilon_r[i]-epsilon_r[i+1]
        # print(K_prime/K)
        cap = pyc.epsilon_0*e_r*(K_prime/K)
        
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
    # print('Co distance',distance_to_ground)
    # print('Co width',width_cap)
    # print('Co vertical',vertical_length_cap)
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
        if np.allclose(k_prime, 1):
            K_prime = ellipkm1(k_prime)
        else:
            K_prime = ellipk(k_prime)
        if i == len(epsilon_r)-1:
            e_r = epsilon_r[i]-1
        else:
            e_r = epsilon_r[i]-epsilon_r[i+1]
        cap = pyc.epsilon_0*(e_r)*(K_prime/K)#**(-1)
        capacitance_contributions.append(cap*vertical_length_cap*2)
    # print(capacitance_contributions)
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
    # print(f' Geometric: {geometric*1e9} nH/m')
    # print(f' Kinetic: {kinetic*1e9*length_ind} nH/m')
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
    # print(f' Ltot: {Lg*1e9} nH')
    # print(f' Cg: {1e15*Cg} fF')
    # print(f' Cc: {1e15*Cc} fF')
    return 1/(2*np.pi*np.sqrt(Lg*(Cc+Cg)))


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
