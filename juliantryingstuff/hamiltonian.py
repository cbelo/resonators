
import numpy as np
import scipy.linalg as sp

def make_hamiltonian(N, wr, Ji, Jx, Jprime, SNN = True):
    H = np.zeros((N,N))

    for i in range(N): 
        H[i,i] = wr

    for j in range(int(N/2)): 
        H[2*j, 2*j+1] = Ji
        H[2*j +1, 2*j] = Ji

    for m in range(int(N/2 -1)):
        H[2*m+2, 2*m +1] = Jx
        H[2*m+1, 2*m+2] = Jx

    if SNN: 
        for x in range(N-2): 
            H[x, x+2] = Jprime
            H[x+2, x] = Jprime
    return H


def make_hamiltonian_diffSNN(N, wr, Ji, Jx, Jprime, SNN = True):
    H = np.zeros((N,N))

    for i in range(N): 
        H[i,i] = wr

    for j in range(int(N/2)): 
        H[2*j, 2*j+1] = Ji
        H[2*j +1, 2*j] = Ji

    for m in range(int(N/2 -1)):
        H[2*m+2, 2*m +1] = Jx
        H[2*m+1, 2*m+2] = Jx

    if SNN: 
        for x in range(1,N-3): 
            H[x, x+2] = Jprime
            H[x+2, x] = Jprime
    return H


def make_hamiltonian_Cinput(N, C0, Ci, Cx, Cprime, L0, SNN = True):
    Csum = Ci + Cx + C0 + 2*Cprime
    w0 = 1/(np.sqrt(Csum*L0))
    Ji = w0/2 * Ci/Csum
    Jx = w0/2 * Cx/Csum
    Jprime = w0/2 * Cprime/Csum
    H = np.zeros((N,N))

    for i in range(N): 
        H[i,i] = w0

    for j in range(int(N/2)): 
        H[2*j, 2*j+1] = Ji
        H[2*j +1, 2*j] = Ji

    for m in range(int(N/2 -1)):
        H[2*m+2, 2*m +1] = Jx
        H[2*m+1, 2*m+2] = Jx

    if SNN: 
        for x in range(N-2): 
            H[x, x+2] = Jprime
            H[x+2, x] = Jprime
    return H

def add_TNN(N, H, J_TNN): 
    for i in range(N-3): 
        H[i, i+3] = J_TNN
        H[i+3,i] = J_TNN

    return H


def make_hamiltonian_withLix(N, C0, Ci, Cx, Cprime, L0, Li, Lx, SNN = True):
    Csum = Ci + Cx + C0 + 2*Cprime
    Lsum = 1/(1/L0+1/Li+1/Lx)
    w0 = 1/(np.sqrt(Csum*Lsum))
    # w0 = 1/(np.sqrt(Csum*L0))
    Ji = w0/2 * Ci/Csum
    Jx = w0/2 * Cx/Csum
    Gi = w0/2 * Lsum/Li
    Gx = w0/2 * Lsum/Lx
    Jprime = w0/2 * Cprime/Csum
    H = np.zeros((N,N))

    for i in range(N): 
        H[i,i] = w0

    for j in range(int(N/2)): 
        H[2*j, 2*j+1] = Ji - Gi
        H[2*j +1, 2*j] = Ji - Gi

    for m in range(int(N/2 -1)):
        H[2*m+2, 2*m +1] = Jx - Gx
        H[2*m+1, 2*m+2] = Jx - Gx

    if SNN: 
        for x in range(N-2): 
            H[x, x+2] = Jprime
            H[x+2, x] = Jprime
    return H

def make_hamiltonian_withLixLprime(N, C0, Ci, Cx, Cprime, L0, Li, Lx, Lprime, SNN = True):
    Csum = Ci + Cx + C0 + 2*Cprime
    Lsum = 1/(1/L0+1/Li+1/Lx)
    # w0 = 1/(np.sqrt(Csum*Lsum))
    w0 = 1/(np.sqrt(Csum*L0))
    Ji = w0/2 * Ci/Csum
    Jx = w0/2 * Cx/Csum
    Jprime = w0/2 * Cprime/Csum
    Gi = w0/2 * Lsum/Li
    Gx = w0/2 * Lsum/Lx
    Gprime = w0/2 * Lsum/Lprime
    H = np.zeros((N,N))

    for i in range(N): 
        H[i,i] = w0

    for j in range(int(N/2)): 
        H[2*j, 2*j+1] = Ji - Gi
        H[2*j +1, 2*j] = Ji - Gi

    for m in range(int(N/2 -1)):
        H[2*m+2, 2*m +1] = Jx - Gx
        H[2*m+1, 2*m+2] = Jx - Gx

    if SNN: 
        for x in range(N-2): 
            H[x, x+2] = Jprime - Gprime
            H[x+2, x] = Jprime - Gprime
    return H






def make_hamiltonian_withLix_resinput(N, f, C0, Ci, Cx, Cprime, L0, Li, Lx, SNN = True):
    Csum = Ci + Cx + C0 + 2*Cprime
    Lsum = 1/(1/L0+1/Li+1/Lx)
    # w0 = 1/(np.sqrt(Csum*Lsum))
    w0 = np.pi*2*f
    Ji = w0/2 * Ci/Csum
    Jx = w0/2 * Cx/Csum
    Gi = w0/2 * Li/Lsum
    Gx = w0/2 * Lx/Lsum
    Jprime = w0/2 * Cprime/Csum
    H = np.zeros((N,N))

    for i in range(N): 
        H[i,i] = w0

    for j in range(int(N/2)): 
        H[2*j, 2*j+1] = Ji - Gi
        H[2*j +1, 2*j] = Ji - Gi

    for m in range(int(N/2 -1)):
        H[2*m+2, 2*m +1] = Jx - Gx
        H[2*m+1, 2*m+2] = Jx - Gx

    if SNN: 
        for x in range(N-2): 
            H[x, x+2] = Jprime
            H[x+2, x] = Jprime
    return H


def build_omega(C0, L0, Cw):
    '''
    Gives the resonance frequency of the resonators in the SSH chain.
    '''
    omega = np.sqrt(1/((L0**(-1))**(-1)*(C0 + Cw)))
    return omega

def C_SNN_contribution(index, C_SNN, N):
    '''
    Returns the contribution of the second neighbour coupling to the index-th resonator.
    '''
    Csnn = 0
    if index > 0:
        Csnn += C_SNN[index-1]
    if index < 2*N-1:
        Csnn += C_SNN[index+1]
    return Csnn


def make_chain(N, L0, C0, Cw, Cv, only_eigval = False, C_SNN = None, diag_elems = False):
    '''
    Returns the eigenvalues and eigenvectors of the SSH chain with Lv and Cw.

    N: number of unit cells (resonators/2)
    '''
    if (N +1) != len(Cw):
        raise ValueError('The length of Cw must be equal to N+1.')
    if N != len(Cv):
        raise ValueError('The length of Cv must be equal to N.')
    if 2*N != len(L0):
        raise ValueError('The length of L0 must be equal to 2*N.')
    if 2*N != len(C0):
        raise ValueError('The length of C0 must be equal to 2*N.')
    if (2*N) != len(C_SNN):
        raise ValueError('The length of C_SNN must be equal to 2*N.')
    if C_SNN is None:
        C_SNN = np.zeros(2*N-2)
    

    #Build the Hamiltonian
    ham = np.zeros((2*N, 2*N))
    #On-site energies
    for i in range(2*N):
        Csnn = C_SNN_contribution(i, C_SNN, N)

        if i%2 == 0:
            omega = build_omega(C0[i], L0[i], (Cw[i//2] + Cv[i//2]+Csnn))
            # print('freq = ', omega/(2*np.pi*1e9), 'GHz')

        else:
            omega = build_omega(C0[i], L0[i], (Cw[(i-1)//2 +1]+Cv[(i-1)//2]+Csnn))
        
        ham[i, i] = omega

    if diag_elems:
        print(np.diag(ham)/(2*np.pi*1e9))
    
    #Coupling terms
    for i in range(2*N-1):
        Csnn_i = C_SNN_contribution(i, C_SNN, N)
        Csnn_i1 = C_SNN_contribution(i+1, C_SNN, N)
        if i%2 == 0:
            ham[i, i+1] = np.sqrt(build_omega(C0[i], L0[i], (Cw[i//2]+Cv[i//2]+Csnn_i))*build_omega(C0[i+1], L0[i+1], (Cw[i//2+1]+Cv[i//2]+Csnn_i1)))/2 * (Cv[i//2]/np.sqrt((C0[i]+Cw[i//2]+Cv[i//2]+Csnn_i)*(C0[i+1]+Cw[i//2+1]+Cv[i//2]+Csnn_i1)))
            ham[i+1,i] =  np.sqrt(build_omega(C0[i], L0[i], (Cw[i//2]+Cv[i//2]+Csnn_i))*build_omega(C0[i+1], L0[i+1], (Cw[i//2+1]+Cv[i//2]+Csnn_i1)))/2 * (Cv[i//2]/np.sqrt((C0[i]+Cw[i//2]+Cv[i//2]+Csnn_i)*(C0[i+1]+Cw[i//2+1]+Cv[i//2]+Csnn_i1)))
        else:
            ham[i, i+1] = np.sqrt(build_omega(C0[i], L0[i], (Cw[(i+1)//2]+Cv[(i+1)//2-1]+Csnn_i))*build_omega(C0[i+1], L0[i+1], (Cw[(i+1)//2]+Cv[(i+1)//2]+Csnn_i1) ))/2 * Cw[(i+1)//2]/np.sqrt((C0[i]+Cw[(i+1)//2]+Cv[(i+1)//2-1]+Csnn_i)*(C0[i+1]+Cw[(i+1)//2]+Cv[(i+1)//2]+Csnn_i1))
            ham[i+1,i] =  np.sqrt(build_omega(C0[i], L0[i], (Cw[(i+1)//2]+Cv[(i+1)//2-1]+Csnn_i))*build_omega(C0[i+1], L0[i+1], (Cw[(i+1)//2]+Cv[(i+1)//2]+Csnn_i1)))/2 * Cw[(i+1)//2]/np.sqrt((C0[i]+Cw[(i+1)//2]+Cv[(i+1)//2-1]+Csnn_i)*(C0[i+1]+Cw[(i+1)//2]+Cv[(i+1)//2]+Csnn_i1))

    #Second neighbour coupling
    if True: #np.allclose(C_SNN, 0):
        for i in range(2, 2*N):
            j = i-2
            Csnn_i = C_SNN_contribution(i, C_SNN, N)
            Csnn_j = C_SNN_contribution(j, C_SNN, N)
            if i%2 == 0:
                wi = build_omega(C0[i], L0[i], (Cw[i//2]+Cv[i//2]+Csnn_i))
                wj = build_omega(C0[j], L0[j], (Cw[j//2]+Cv[j//2]+Csnn_j))
                omega = np.sqrt(wi*wj)
                denominator = 2*np.sqrt((C0[i]+Cw[i//2]+Cv[i//2]+Csnn_i)*(C0[j]+Cw[j//2]+Cv[j//2]+Csnn_j))
                ham[i, j] = C_SNN[j]*omega/denominator
                ham[j, i] = C_SNN[j]*omega/denominator
            else:
                wi = build_omega(C0[i], L0[i], (Cw[(i+1)//2]+Cv[(i+1)//2-1]+Csnn_i))
                wj = build_omega(C0[j], L0[j], (Cw[(j+1)//2]+Cv[(j+1)//2-1]+Csnn_j))
                omega = np.sqrt(wi*wj)
                denominator = 2*np.sqrt((C0[i]+Cw[(i+1)//2]+Cv[(i+1)//2-1]+Csnn_i)*(C0[j]+Cw[(j+1)//2]+Cv[(j+1)//2-1]+Csnn_j))
                ham[i, j] = C_SNN[j]*omega/denominator
                ham[j, i] = C_SNN[j]*omega/denominator



    # np.set_printoptions(precision=2, linewidth = 140)
    # print(ham)
    #Get the eigenvalues and eigenvectors
    if only_eigval:
        eigval = sp.eigvalsh(ham)
        # print(eigval/(2*np.pi*1e9))
        return eigval
    else:
        eigval, eigvec = sp.eig(ham)
        return eigval, eigvec.T, ham #Transpose to get the eigenvectors in rows so eigvec[i] is the eigenvector of eigval[i]
