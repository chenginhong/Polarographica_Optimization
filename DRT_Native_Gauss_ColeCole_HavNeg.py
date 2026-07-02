# -*- coding: utf-8 -*-
"""
#=======================================================================================================================
#=======================================================================================================================
#Author:       Tim Tichter
#Year:         2022
#Function:     DRT functions of POLAROGRAPHICA
#=======================================================================================================================
#=======================================================================================================================
#     Here, three additional ways of computing the DRT function are defined
#     The work behind this function was carried out during Tim Tichters stay at DTU
#     Technical University of Denmark in 2022.
#     The spike-DRT is the simplest version of computing a DRT (the basis function is a
#     Dirac Deltafunction). The other two types of DRTs use either a Gaussian distribution or
#     the "true" analytical DRT representation for the Cole-Cole constant phase elemnt.
#=======================================================================================================================
#=======================================================================================================================
"""

import numpy                           as np
from scipy.optimize                    import nnls
from math                              import ceil,floor

##############################################################################################
#  The following two functions are used for the interconversion of the natural and the 
#  decadic logarithm.
##############################################################################################
def ln_to_log(x):
    return x/np.log(10)

def log_to_ln(x):
    return x*np.log(10)

##############################################################################################
#  The following function is the "two-parameter-four-quadrant-inverse-tangent".
#  It is used to unabmiguously define angles computed by arctan, i.e. it defines points
#  also in the left half-plane of the complex plane. For this purpose, it adds \pi, if
#  the "real" argument of the arctan becomes negative. This corresponds to negative real parts.
#  Note that the sign of the imaginary part was removed in the derivation of the actual 
#  DRT-formula by setting arctan(-x) = -arctan(x). If a negative argument occurs, it is 
#  IN TIS CASE from negative real parts (cf. derivation).
##############################################################################################
def ARCTAN(INPUT):
    RESULT = []
    for i in range(len(INPUT)):
        if INPUT[i] >= 0:
            RESULT.append(np.arctan(INPUT[i]))
        else:
            RESULT.append(np.arctan(INPUT[i]) + np.pi)
    return np.array(RESULT)

##############################################################################################
#  In this step, the Gaussian and Cole-Cole distribution matrices are generated.
#  This is required to put a smoothing into the final fitting step. The Gaussian and
#  Cole-Cole DRTs are normalized - i.e. their integral is equal to one.
##############################################################################################

def build_Gauss(LOGtau, FWHM):
    c = FWHM/2.35482
    X_Gauss = np.zeros((len(LOGtau),len(LOGtau)))
    for m in range(len(LOGtau)):
        X_Gauss[::,m] = (np.exp(-((LOGtau[m]-LOGtau)**2)/(2*c**2)))/(2*np.pi*c**2)**0.5
    return X_Gauss


def build_ColeCole(LOGtau, GAMMA):
    Taus = 2.7182818**LOGtau
    X_ColeCole = np.zeros((len(LOGtau),len(LOGtau)))
    for m in range(len(LOGtau)):
        X_ColeCole[::,m] = (1/(2*np.pi))*np.sin(np.pi*(1-GAMMA))/(np.cosh((GAMMA)*np.log(Taus/Taus[m]))  -  np.cos(np.pi*(1-GAMMA))   )
    return X_ColeCole


def build_HavriliakNegami(LOGtau, ALPHA, BETA):
    Taus = 2.7182818**LOGtau
    X_HavriliakNegami = np.zeros((len(LOGtau),len(LOGtau)))
    for m in range(len(LOGtau)):
        THETA = ARCTAN(INPUT = np.sin(np.pi*ALPHA)/( (Taus/Taus[m])**ALPHA  +  np.cos(ALPHA*np.pi)))
        X_HavriliakNegami[::,m] = (1/(np.pi))*( (((Taus/Taus[m])**(ALPHA*BETA))*np.sin(BETA*THETA)) / 
                                ( 1 + ((Taus/Taus[m])**(2*ALPHA)) + 2*((Taus/Taus[m])**(ALPHA))*np.cos(np.pi*ALPHA)   )**(BETA/2) )
    return X_HavriliakNegami
  

  
    
################################################################################################
#  Here, a set of time constants is generated on the base of the frequencied which are provided
#  with the experimental dataset.
################################################################################################    
    
def make_tau(freq_raw, low_ext, high_ext, resol):
    taumax      = ceil(np.max(np.log10(1./freq_raw))) + high_ext
    taumin      = floor(np.min(np.log10(1./freq_raw)))- np.abs(low_ext)
    tau         = np.logspace(taumin,taumax,resol*len(freq_raw))
    return tau    


################################################################################################
#  Here, the RC Kernel-Matrix is generated. This matrix is the key for computing the DRT later on
################################################################################################    

def build_RC(freq_raw, tau):
    X_RC   = np.zeros((len(freq_raw),len(tau)),dtype=np.complex128)
    dlntau = np.log(tau[1]/tau[0])
    for m in range(len(freq_raw)):
        X_RC[m,::] = dlntau/(1+2j*np.pi*freq_raw[m]*tau)
    return X_RC


################################################################################################
# The following function will generate a template for the regularization matrix. This is
# important, since reconstructing the DRT is an inherently ill-posed problem. The penalty
# term, used for the regularization is linear in the most simple case. In this scenario,
# each entry in the DRT function will introduce an equally weighted cost in the fit.
# Optinally, the penalty can be set to a first or second order derivative of the fitting
# result. In case of a first-order derivative, abrupt changes (i.e. slopes) in the result
# will introduce a stronger cost. In case of a second order derivative, vertex-points will
# itroduce a stronger penalty. In any case, the regularization matris is subsequently multiplied
# with the regularization parameter. 
################################################################################################ 
    
def build_RegMat_Template(Size, delta, Reg_Type):
    if Reg_Type     == '1st-order':
        Reg_Template  = np.zeros((Size, Size))
        for i in range(Size):
            if i > 0 and i < (Size-1):
                Reg_Template[i, i+1] =  1/(2*delta)
                Reg_Template[i, i-1] = -1/(2*delta)
            if i == 0:
                Reg_Template[0, 0]   = -3/(2*delta)
                Reg_Template[0, 1]   =  4/(2*delta)
                Reg_Template[0, 2]   = -1/(2*delta)
            if i == (len(Reg_Template[::,0])-1):
                Reg_Template[-1,-3] =  1/(2*delta)
                Reg_Template[-1,-2] = -4/(2*delta)
                Reg_Template[-1,-1] =  3/(2*delta)
    elif Reg_Type   == '2nd-order':
        Reg_Template    = np.zeros((Size, Size))
        for i in range(Size):
            if i > 0 and i < (Size-1):
                Reg_Template[i, i+1] =  1/(delta**2)
                Reg_Template[i, i  ] = -2/(delta**2)
                Reg_Template[i, i-1] =  1/(delta**2)
            if i == 0:
                Reg_Template[0, 0]  =  2/(delta**2)
                Reg_Template[0, 1]  = -5/(delta**2)
                Reg_Template[0, 2]  =  4/(delta**2)
                Reg_Template[0, 3]  = -1/(delta**2)
            if i == (Size-1):
                Reg_Template[-1,-4] = -1/(delta**2)
                Reg_Template[-1,-3] =  4/(delta**2)
                Reg_Template[-1,-2] = -5/(delta**2)
                Reg_Template[-1,-1] =  2/(delta**2)
    elif Reg_Type   == 'linear':
        Reg_Template    = np.identity(Size)
    return Reg_Template


################################################################################################
#  This function computes the native DRT. This requires regularization. We use Tikhonov-regularization
#  for this purpose. The resulting Matrix-vector equation gives a set of linear equations
#  which are solved by the Lawson-Hansen NNLS (non-negative least-square) solver.
#  Regularization is achieved by adding a linear penalty term (so no derivative penalty).
#  Essentially, we do
#
#  | X_RC | |DRT|  = |Z|
#  | e*U  |          |0|
#  where U is the unity matrix, e the regularization parameter and 0 he zero vector.
#  For the Real part, we take |X_RC.real| and Z.real
#  For the Imaginary part, we take |X_RC.imag| and Z.imag
#  and for the combi-fit, we do
#  | X_RC_re | |DRT|  = |Z_re|
#  | X_RC_im |          |Z_im| 
#  | e*U     |          |0|
################################################################################################ 

def DRT_Native(freq_raw, Z_raw, damp, low_ext, high_ext, resol, DRT_Type, Reg_Type):
    tau             = make_tau(freq_raw, low_ext, high_ext, resol)
    dlntau          = np.log(tau[1]/tau[0])
    X_RC_classic    = build_RC(freq_raw, tau)
    X_RC            = X_RC_classic
    RegulMat        = build_RegMat_Template(X_RC.shape[1], dlntau, Reg_Type)
    if DRT_Type == "real":
        Kernel_Matrix = np.vstack((X_RC.real,damp*RegulMat))
        gamma_raw     = nnls(Kernel_Matrix, np.concatenate([Z_raw.real,np.zeros(len(tau))]))[0]
    if DRT_Type == "imag":
        Kernel_Matrix = np.vstack((X_RC.imag,damp*RegulMat))
        gamma_raw     = nnls(Kernel_Matrix, np.concatenate([Z_raw.imag,np.zeros(len(tau))]))[0]
    if DRT_Type == "comb":
        Kernel_Matrix = np.vstack((np.vstack((X_RC.real,X_RC.imag)),damp*RegulMat)) 
        gamma_raw     = nnls(Kernel_Matrix, np.concatenate([np.concatenate([Z_raw.real,Z_raw.imag]),np.zeros(len(tau))]))[0]
    return np.log10(tau), gamma_raw, Kernel_Matrix, X_RC_classic, X_RC, dlntau


################################################################################################
#  This function computes the Gaussian DRT. This again requires regularization. We use Tikhonov
#  -regularization for this purpose. The resulting Matrix-vector equation gives a set of linear 
#  equations which are solved by the Lawson-Hansen NNLS (non-negative least-square) solver.
#  Regularization is achieved by adding a linear penalty term (so no derivative penalty).
#  Essentially, we take
#  X_RC_classic  --> The kernel matrix from the classical DRT
#  X_Gauss       --> A Matrix which ntroduces a broadening of each DRT signal according to
#                    a normalized Gauss function (its integral is one).
#  X_RC =  dot(X_RC_classic, X_Gauss) --> The dot product of both matrices. Subsequently, we do
#  ----------------------------------------------------------
#  | X_RC | |DRT|  = |Z|
#  | e*U  |          |0|
#  ----------------------------------------------------------
#  where U is the unity matrix, e the regularization parameter and 0 he zero vector.
#  For the Real part, we take |X_RC.real| and Z.real
#  For the Imaginary part, we take |X_RC.imag| and Z.imag
#  and for the combi-fit, we do
#  | X_RC_re | |DRT|  = |Z_re|
#  | X_RC_im |          |Z_im| 
#  | e*U     |          |0|
################################################################################################ 

def DRT_Gaussian_RBF(freq_raw, Z_raw, damp, low_ext, high_ext, resol, DRT_Type, FWHM, Reg_Type):
    tau             = make_tau(freq_raw, low_ext, high_ext, resol)
    dlntau          = np.log(tau[1]/tau[0])
    X_RC_classic    = build_RC(freq_raw, tau)
    X_Gauss         = build_Gauss(np.log(tau), np.log(10)*FWHM)
    X_RC            = np.dot(X_RC_classic, X_Gauss)
    RegulMat        = build_RegMat_Template(X_RC.shape[1], dlntau, Reg_Type)
    if DRT_Type == "real":
        Kernel_Matrix = np.vstack((X_RC.real,damp*RegulMat))
        gamma_raw     = nnls(Kernel_Matrix, np.concatenate([Z_raw.real,np.zeros(len(tau))]))[0]
    if DRT_Type == "imag":
        Kernel_Matrix = np.vstack((X_RC.imag,damp*RegulMat))
        gamma_raw     = nnls(Kernel_Matrix, np.concatenate([Z_raw.imag,np.zeros(len(tau))]))[0]
    if DRT_Type == "comb":
        Kernel_Matrix = np.vstack((np.vstack((X_RC.real,X_RC.imag)),damp*RegulMat)) 
        gamma_raw     = nnls(Kernel_Matrix, np.concatenate([np.concatenate([Z_raw.real,Z_raw.imag]),np.zeros(len(tau))]))[0]
    gamma_fine = np.dot(X_Gauss, gamma_raw)
    return np.log10(tau), gamma_raw, gamma_fine, Kernel_Matrix, X_Gauss, X_RC_classic, X_RC, dlntau


###################################################################################################
#  Here, essentially the same is done as for the Gaussian DRT. The difference is that the analytical
#  DRT expression from Cole and Cole [1] is used as radial basis function. 
#  NOTE THAT THIS EXPRESSION IS VALID FOR A CPE OF THE FORM
#
#  Z(\omega)-Z(\infty)/(Z(\omega)-Z(\infty)) = 1/(1 + (i\omega\tau_{0})^{1-\phi} )
#  
#  where \phi is the ideality-parameter of the CPE. We have taken (1-\phi) = \gamma (not to confuse
#  with the DRT-\gamma(\tau))). In this manner, we see that if \gamma = 1, we have an ideal capacitor.
#  and 0 < \gamma <= 1.
#  
#  [1] K. S. Cole, R. H. Cole, Dispersion and Absorption in Dielectrics I. Alternating Current 
#      Characteristics, The Journal of Chemical Physics 9 (1941)295 341{351. doi:10.1063/1.1750906.
####################################################################################################

def DRT_ColeCole_RBF(freq_raw, Z_raw, damp, low_ext, high_ext, resol, DRT_Type, GAMMA_Cole, Reg_Type):
    tau             = make_tau(freq_raw, low_ext, high_ext, resol)
    dlntau          = np.log(tau[1]/tau[0])
    X_RC_classic    = build_RC(freq_raw, tau)
    X_ColeCole      = build_ColeCole(np.log(tau), GAMMA_Cole)   # Natural-log of the Taus
    X_RC            = np.dot(X_RC_classic, X_ColeCole)
    RegulMat        = build_RegMat_Template(X_RC.shape[1], dlntau, Reg_Type)
    if DRT_Type == "real":
        Kernel_Matrix = np.vstack((X_RC.real,damp*RegulMat))
        gamma_raw     = nnls(Kernel_Matrix, np.concatenate([Z_raw.real,np.zeros(len(tau))]))[0]
    if DRT_Type == "imag":
        Kernel_Matrix = np.vstack((X_RC.imag,damp*RegulMat))
        gamma_raw     = nnls(Kernel_Matrix, np.concatenate([Z_raw.imag,np.zeros(len(tau))]))[0]
    if DRT_Type == "comb":
        Kernel_Matrix = np.vstack((np.vstack((X_RC.real,X_RC.imag)),damp*RegulMat)) 
        gamma_raw     = nnls(Kernel_Matrix, np.concatenate([np.concatenate([Z_raw.real,Z_raw.imag]),np.zeros(len(tau))]))[0]
    gamma_fine = np.dot(X_ColeCole, gamma_raw)
    return np.log10(tau), gamma_raw, gamma_fine, Kernel_Matrix, X_ColeCole, X_RC_classic, X_RC, dlntau



###################################################################################################
#  Here, the Havriliak-Negami DRT is defined. For this purpose, the analytical
#  DRT expression of the Havriliak Negami relaxation is used as radial basis function. 
#  NOTE THAT THIS EXPRESSION IS VALID FOR A CPE OF THE FORM
#
#  Z(\omega)-Z(\infty)/(Z(\omega)-Z(\infty)) = 1/((1 + (i\omega\tau_{0})^{alpha} )^beta)
#  -------------------------------------------------------------------------------------------------
#  where alpha and beta tune the width and skewness of the relaxation function.
#  in case of     beta = 1, alpha = 1   --->   Classical spike DRT (Debye relaxation)
#  in case of     beta = 1, alpha < 1   --->   Cole-Cole relaxation
#  in case of     beta < 1, alpha = 1   --->   Cole-Davidson relaxation (geed for i.e. RDE)
#  in case of     beta < 1, alpha < 1   --->   Actual Havriliak-Negami model
#  -------------------------------------------------------------------------------------------------
####################################################################################################

def DRT_HavriliakNegami_RBF(freq_raw, Z_raw, damp, low_ext, high_ext, resol, DRT_Type, ALPHA_HAVNEG, BETA_HAVNEG, Reg_Type):
    tau             = make_tau(freq_raw, low_ext, high_ext, resol)
    dlntau          = np.log(tau[1]/tau[0])
    X_RC_classic    = build_RC(freq_raw, tau)
    X_HavNeg        = build_HavriliakNegami(np.log(tau), ALPHA_HAVNEG, BETA_HAVNEG)    # Natural-log of the Taus
    X_RC            = np.dot(X_RC_classic, X_HavNeg)
    RegulMat        = build_RegMat_Template(X_RC.shape[1], dlntau, Reg_Type)
    if DRT_Type == "real":
        Kernel_Matrix = np.vstack((X_RC.real,damp*RegulMat))
        gamma_raw     = nnls(Kernel_Matrix, np.concatenate([Z_raw.real,np.zeros(len(tau))]))[0]
    if DRT_Type == "imag":
        Kernel_Matrix = np.vstack((X_RC.imag,damp*RegulMat))
        gamma_raw     = nnls(Kernel_Matrix, np.concatenate([Z_raw.imag,np.zeros(len(tau))]))[0]
    if DRT_Type == "comb":
        Kernel_Matrix = np.vstack((np.vstack((X_RC.real,X_RC.imag)),damp*RegulMat)) 
        gamma_raw     = nnls(Kernel_Matrix, np.concatenate([np.concatenate([Z_raw.real,Z_raw.imag]),np.zeros(len(tau))]))[0]
    gamma_fine = np.dot(X_HavNeg, gamma_raw)
    return np.log10(tau), gamma_raw, gamma_fine, Kernel_Matrix, X_HavNeg, X_RC_classic, X_RC, dlntau

###################################################################################################
# The following function will perform a Kramers-Kronig like back-multiplication of a given
# RC-matrix and a DRT. This means, Real-part of the matrix X_RC will be multiplied with the
# DRT computed from the imaginary part of the impedance and the imaginary part of X_RC will
# be multiplied with the DRT computed from the real part of the impedance to give the imaginary
# part of the back-computed impedance. This can be seen as a Kramers-Kronig-check. Thus,
# if the back-computed Impedance data matches with the experimental data, the data fulfils the
# criteria of linearity, causality and time-invarianc --> i.e. is good data :)
###################################################################################################

def KramersKronig_Multiplication(MATRIX, Real_DRT, Imag_DRT):
    Z_real_back = np.dot(MATRIX.real, Imag_DRT)
    Z_imag_back = np.dot(MATRIX.imag, Real_DRT)
    Z_KK_back   = Z_real_back + 1j*Z_imag_back
    return Z_KK_back

def AntiKramersKronig_Multiplication(MATRIX, Real_DRT, Imag_DRT):
    Z_real_back     = np.dot(MATRIX.real, Real_DRT)
    Z_imag_back     = np.dot(MATRIX.imag, Imag_DRT)
    Z_antiKK_back   = Z_real_back + 1j*Z_imag_back
    return Z_antiKK_back

def CombiKramersKronig_Multiplication(MATRIX, Comb_DRT):
    Z_real_back     = np.dot(MATRIX.real, Comb_DRT)
    Z_imag_back     = np.dot(MATRIX.imag, Comb_DRT)
    Z_CombiKK_back  = Z_real_back + 1j*Z_imag_back
    return Z_CombiKK_back