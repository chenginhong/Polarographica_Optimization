# -*- coding: utf-8 -*-
"""
#=======================================================================================================================
#=======================================================================================================================
#Authosr:      Tim Tichter, Jonathan Schneider
#Date:         2019.07.16
#Function:     Pythonic translation of DRT-tools which uses the Lawson-Hanson Algorithm
#              instead of the original QuadProg Algorithm for fitting
#=======================================================================================================================
#=======================================================================================================================
IMPORTANT:     This function is a pythonic version/translation of the great open-source software DRTtools 
               (originally written in Matlab). The main difference is, that we use the  Lawson-Hanson (NNLS) algorithm 
               for optimization instead of quadprog. WHEN USING THIS FUNCTION, CHECK THE FOLLOWING ORIGINAL WORK BY!
               T. H. Wan, M. Saccoccio, C. Chen and F. Ciucci DOI: 10.1016/j.electacta.2015.09.097 as well as the 
               following web-pages https://ciucci.org/project/drt/\nhttps://github.com/ciuccislab/DRTtools")
#=======================================================================================================================
#=======================================================================================================================
"""

import numpy                           as np
from scipy.optimize                    import nnls
from math                              import ceil,floor
from scipy.integrate                   import quad
from scipy.linalg                      import toeplitz


#============================================================================================================   
#                 -------> DRT_TOOLS_NNLS_DRT <-----------
#============================================================================================================ 

#============================================================================================================   
#    -------> define all DRT Tools functions in python 2.7 language except quadprog <-----------
#       as there is not quadprog in python this transfom uses nnls algorthm in the end
#============================================================================================================ 

        
def quad_format_combined(A_re,A_im,b_re,b_im,M_re,M_im,damp): 
    H = 2*(0.5*(np.dot(A_re.T,A_re)+np.dot(A_im.T,A_im))+damp*M_re)
    c = -2*0.5*(np.dot(b_im.T,A_im)+np.dot(b_re.T,A_re))
    return H , c

def quad_format(A,b,M,damp):
    H = 2*(np.dot(A.T,A)+damp*M)
    c = -2*np.dot(b.T,A)
    return H , c    

##############################################################################
#map array to gamma
##############################################################################

def map_array_to_gamma(freq_map, freq_coll, x, epsilon, rbf_type):
    if rbf_type == "gaussian":
        rbf = lambda y,y0: np.exp(-(epsilon*(y-y0))**2)
    elif rbf_type == "C0_matern":
        rbf = lambda y,y0: np.exp(-abs(epsilon*(y-y0)))
    elif rbf_type == "C2_matern":
        rbf = lambda y,y0: np.exp(-abs(epsilon*(y-y0)))*(1+abs(epsilon*(y-y0)))
    elif rbf_type == "C4_matern":
        rbf = lambda y,y0: 1/3.0*np.exp(-abs(epsilon*(y-y0)))*(3+3*abs(epsilon*(y-y0))+abs(epsilon*(y-y0))**2)
    elif rbf_type == "C6_matern":
        rbf = lambda y,y0: 1/15.0*np.exp(-abs(epsilon*(y-y0)))*(15+15*abs(epsilon*(y-y0))+6*abs(epsilon*(y-y0))**2+abs(epsilon*(y-y0))**3)          
    elif rbf_type == "inverse_quadratic":
        rbf = lambda y,y0: 1.0/(1+(epsilon*(y-y0))**2)
    else:
        print('ERROR - Unexpected RBF input at map_array_to_gamma')
    y0  = -np.log(freq_coll)
    out_gamma = np.zeros(len(freq_map))   
    
    
    for iter_freq_map in range(len(freq_map)):
        freq_map_loc = freq_map[iter_freq_map]
        y = -np.log(freq_map_loc)
        rbf_temp = rbf(y,y0)
        out_gamma[iter_freq_map] = np.dot(x.T,rbf_temp)
    
    return out_gamma

##############################################################################
#inner products
##############################################################################
def inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type):
    a = epsilon*np.log(freq_n/freq_m)
    if rbf_type == "gaussian":
        out_IP = epsilon**3*(3-6*a**2+a**4)*np.exp(-(a**2/2))*(np.pi/2)**0.5
        return out_IP
    elif rbf_type == "C0_matern":
        out_IP = epsilon**3*(1+abs(a))*np.exp(-abs(a))
        return out_IP
    elif rbf_type == "C2_matern":
        out_IP = epsilon**3/6.0*(3+3*abs(a)-6*abs(a)**2+abs(a)**3)*np.exp(-abs(a))
        return out_IP
    elif rbf_type == "C4_matern":
        out_IP = epsilon**3/30.0*(45+45*abs(a)-15*abs(a)**3-5*abs(a)**4+abs(a)**5)*np.exp(-abs(a))                     
        return out_IP
    elif rbf_type == "C6_matern":
        out_IP = epsilon**3/140.0*(2835+2835*abs(a)+630*abs(a)**2-315*abs(a)**3-210*abs(a)**4-42*abs(a)**5+abs(a)**7)*np.exp(-abs(a))         
        return out_IP
    elif rbf_type == "inverse_quadratic":
        out_IP = 48*(16+5*a**2*(-8 + a**2))*np.pi*epsilon**3/((4 + a**2)**5)
        return out_IP
    else:
        print('ERROR - Unexpected RBF input at inner_prod_rbf_2')

def inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type):
    a = epsilon*np.log(freq_n/freq_m)

    if rbf_type == "gaussian":
        out_IP = -epsilon*(-1+a**2)*np.exp(-(a**2/2))*(np.pi/2)**0.5
        return out_IP
    elif rbf_type == "C0_matern":
        out_IP = epsilon*(1-abs(a))*np.exp(-abs(a))
        return out_IP
    elif rbf_type == "C2_matern":
        out_IP = epsilon/6.0*(3+3*abs(a)-abs(a)**3)*np.exp(-abs(a))
        return out_IP
    elif rbf_type == "C4_matern":
        out_IP = epsilon/30.0*(105+105*abs(a)+30*abs(a)**2-5*abs(a)**3-5*abs(a)**4-abs(a)**5)*np.exp(-abs(a))                     
        return out_IP
    elif rbf_type == "C6_matern":
        out_IP = epsilon/140.0*(10395 +10395*abs(a)+3780*abs(a)**2+315*abs(a)**3-210*abs(a)**4-84*abs(a)**5-14*abs(a)**6-abs(a)**7)*np.exp(-abs(a))         
        return out_IP
    elif rbf_type == "inverse_quadratic":
        out_IP = 4*epsilon*(4-3*a**2)*np.pi/((4+a**2)**3)
        return out_IP
    else:
        print('ERROR - Unexpected RBF input at inner_prod_rbf.')
    
#############################################################################
#g_i
#############################################################################
def g_i(freq_n, freq_m, epsilon, rbf_type):
    alpha = 2*np.pi*freq_n/freq_m

    if rbf_type == "gaussian":
        rbf = lambda x: np.exp(-(epsilon*x)**2)
    elif rbf_type == "C0_matern":
        rbf = lambda x: np.exp(-abs(epsilon*x))
    elif rbf_type == "C2_matern":
        rbf = lambda x: np.exp(-abs(epsilon*x))*(1+abs(epsilon*x))
    elif rbf_type == "C4_matern":
        rbf = lambda x: (1/3.0)*np.exp(-abs(epsilon*x))*(3+3*abs(epsilon*x)+abs(epsilon*x)**2)
    elif rbf_type == "C6_matern":
        rbf = lambda x: (1/15.0)*np.exp(-abs(epsilon*x))*(15+15*abs(epsilon*x)+6*abs(epsilon*x)**2+abs(epsilon*x)**3)          
    elif rbf_type == "inverse_quadratic":
        rbf = lambda x: 1.0/(1+(epsilon*x)**2)
    else:
        print('ERROR - Unexpected RBF input at g_i')
    
    integrand_g_i = lambda x: 1.0/(1+alpha**2 *np.exp(2*x))*rbf(x)
    
    out_val  = quad(integrand_g_i, -100, 100,epsabs=1.0e-06, epsrel=1.0e-06)[0]#,'RelTol',1E-9,'AbsTol',1e-9);       
    return out_val 


#############################################################################
#g_ii
#############################################################################

def g_ii(freq_n, freq_m, epsilon, rbf_type):
    alpha = 2*np.pi*freq_n/freq_m

    if rbf_type == "gaussian":
        rbf = lambda x: np.exp(-(epsilon*x)**2)
    elif rbf_type == "C0_matern":
        rbf = lambda x: np.exp(-abs(epsilon*x))
    elif rbf_type == "C2_matern":
        rbf = lambda x: np.exp(-abs(epsilon*x))*(1+abs(epsilon*x))
    elif rbf_type == "C4_matern":
        rbf = lambda x: (1/3.0)*np.exp(-abs(epsilon*x))*(3+3*abs(epsilon*x)+abs(epsilon*x)**2)
    elif rbf_type == "C6_matern":
        rbf = lambda x: (1/15.0)*np.exp(-abs(epsilon*x))*(15+15*abs(epsilon*x)+6*abs(epsilon*x)**2+abs(epsilon*x)**3)          
    elif rbf_type == "inverse_quadratic":
        rbf = lambda x: 1.0/(1+(epsilon*x)**2)

    else:
        print('ERROR - Unexpected RBF input at g_ii')
    
    integrand_g_ii = lambda x: alpha/(1.0/np.exp(x)+alpha**2 *np.exp(x))*rbf(x)

    out_val  = quad(integrand_g_ii, -100, 100,epsabs=1.0e-06, epsrel=1.0e-06)[0]       
    return out_val 
    
#############################################################################
#compute_L_re
#############################################################################
def compute_L_re(freq):
    tau    = 1.0/freq
    N_freq = len(freq)
    out_L      = np.zeros((N_freq-1,N_freq+2))
    out_L_temp = np.zeros((N_freq-1,N_freq+1))
    
    for p in range(N_freq-1):
        delta_loc          = np.log(tau[p+1]/tau[p])
        out_L_temp[p,p+1]  = -1/delta_loc
        out_L_temp[p,p+2]  =  1/delta_loc
        
    out_L[::,1::] = out_L_temp   
    return out_L

#############################################################################
#compute_L_im
############################################################################
def compute_L_im(freq):
    tau    = 1.0/freq
    N_freq = len(freq)
    
    out_L      = np.zeros((N_freq-1,N_freq+2))
    out_L_temp = np.zeros((N_freq-1,N_freq))
    
    for p in range(N_freq-1):
        delta_loc        = np.log(tau[p+1]/tau[p])
        out_L_temp[p,p]  = -1/delta_loc
        out_L_temp[p,p+1]=  1/delta_loc
        
    out_L[::,2::] = out_L_temp   
    return out_L

#############################################################################
#assemble A
############################################################################

def assemble_A_im(freq, epsilon, rbf_type,L=0):
    std_freq      = np.std(np.diff(np.log(1.0/freq)))
    mean_freq     = np.mean(np.diff(np.log(1.0/freq)))
    R             = np.zeros((1,len(freq)))
    C             = np.zeros((len(freq),1))
    out_A_im_temp = np.zeros(len(freq))
    out_A_im      = np.zeros((len(freq), len(freq)+2))
    if std_freq/mean_freq < 1:  
        for iter_freq_n in range(len(freq)): 
            freq_n     = freq[iter_freq_n]
            freq_m = freq[0]
            C[iter_freq_n, 0] = g_ii(freq_n, freq_m, epsilon, rbf_type) 

        for iter_freq_m in range(len(freq)):
            freq_n = freq[0]
            freq_m = freq[iter_freq_m]
            R[0, iter_freq_m] = g_ii(freq_n, freq_m, epsilon, rbf_type)

        out_A_im_temp = toeplitz(C,R)
    
    else:
        for iter_freq_n in range(len(freq)):
            for iter_freq_m in range(len(freq)):
                freq_n = freq[iter_freq_n] 
                freq_m = freq[iter_freq_m]
                out_A_im_temp[iter_freq_n, iter_freq_m] = g_ii(freq_n, freq_m, epsilon, rbf_type) 
    out_A_im[:, 2::] = out_A_im_temp
    if L==1:
        out_A_im[:,0] = -2*np.pi*(freq[:])         
    return out_A_im


def assemble_A_re(freq, epsilon, rbf_type):
    std_freq      = np.std(np.diff(np.log(1.0/freq)))
    mean_freq     = np.mean(np.diff(np.log(1.0/freq)))
    R             = np.zeros((1,len(freq)))
    C             = np.zeros((len(freq),1))
    out_A_re_temp = np.zeros(len(freq))
    out_A_re      = np.zeros((len(freq), len(freq)+2))
    
    if std_freq/mean_freq < 1:  #(error in frequency difference <1% make sure that the terms are evenly distributed)
        for iter_freq_n in range(len(freq)): 
            freq_n     = freq[iter_freq_n]
            freq_m = freq[0]
            C[iter_freq_n, 0] = g_i(freq_n, freq_m, epsilon, rbf_type) 

        for iter_freq_m in range(len(freq)):
            freq_n = freq[0]
            freq_m = freq[iter_freq_m]
            R[0, iter_freq_m] = g_i(freq_n, freq_m, epsilon, rbf_type)

        out_A_re_temp = toeplitz(C,R)
    
    else:
        for iter_freq_n in range(len(freq)):
            for iter_freq_m in range(len(freq)):
                freq_n = freq[iter_freq_n] 
                freq_m = freq[iter_freq_m]
                out_A_re_temp[iter_freq_n, iter_freq_m] = g_i(freq_n, freq_m, epsilon, rbf_type)
    
    out_A_re[:, 2::] = out_A_re_temp
    out_A_re[:,1] = 1
    return out_A_re
    
##############################################################################
#assemble M
##############################################################################
def assemble_M_re(freq, epsilon, rbf_type, der_used):
    std_freq      = np.std(np.diff(np.log(1.0/freq)))
    mean_freq     = np.mean(np.diff(np.log(1.0/freq)))
    R             = np.zeros((1,len(freq)))
    C             = np.zeros((len(freq),1))
    out_M_re_temp = np.zeros(len(freq))
    out_M_re      = np.zeros((len(freq)+2, len(freq)+2))
    
    if der_used == "1st-order":
        if std_freq/mean_freq < 1:  #%(error in frequency difference <1% make sure that the terms are evenly distributed)
            for iter_freq_n in range(len(freq)):
                freq_n = freq[iter_freq_n]
                freq_m = freq[0]
                C[iter_freq_n, 0] = inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type)
    
            for iter_freq_m in range(len(freq)):
                freq_n = freq[0]
                freq_m = freq[iter_freq_m]
                R[0, iter_freq_m] = inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type)
                
            out_M_re_temp = toeplitz(C,R)

        else:
            for iter_freq_n in range(len(freq)):
                for iter_freq_m in range(len(freq)):
                    freq_n = freq[iter_freq_n]
                    freq_m = freq[iter_freq_m]
                    out_M_re_temp[iter_freq_n, iter_freq_m] = inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type)
                    
    #-------------------------------------------------------------------------------------------------------------------           
       
    if der_used == "2nd-order":
        if std_freq/mean_freq < 1:  #%(error in frequency difference <1  #% make sure that the terms are evenly distributed)
            for iter_freq_n in range(len(freq)):
                freq_n = freq[iter_freq_n]
                freq_m = freq[0]
                C[iter_freq_n, 0] = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type)
                    
            for iter_freq_m in range(len(freq)):

                freq_n = freq[0]
                freq_m = freq[iter_freq_m]
                R[0, iter_freq_m] = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type)
                
            out_M_re_temp = toeplitz(C,R)

        else: #%if log of tau is not evenly distributed
            for iter_freq_n in range(len(freq)):
                for iter_freq_m in range(len(freq)):
                    freq_n = freq[iter_freq_n]
                    freq_m = freq[iter_freq_m]
                    out_M_re_temp[iter_freq_n, iter_freq_m] = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type)

    out_M_re[2::, 2::] = out_M_re_temp
    
    return out_M_re


def assemble_M_im(freq, epsilon, rbf_type, der_used):
    std_freq      = np.std(np.diff(np.log(1.0/freq)))
    mean_freq     = np.mean(np.diff(np.log(1.0/freq)))
    R             = np.zeros((1,len(freq)))
    C             = np.zeros((len(freq),1))
    out_M_im_temp = np.zeros(len(freq))
    out_M_im      = np.zeros((len(freq)+2, len(freq)+2))

    if der_used == "1st-order":
        if std_freq/mean_freq < 1:  #%(error in frequency difference <1% make sure that the terms are evenly distributed)
            for iter_freq_n in range(len(freq)):
                freq_n = freq[iter_freq_n]
                freq_m = freq[0]
                C[iter_freq_n, 0] = inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type)
    
            for iter_freq_m in range(len(freq)):
                freq_n = freq[0]
                freq_m = freq[iter_freq_m]
                R[0, iter_freq_m] = inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type)
            
            out_M_im_temp = toeplitz(C,R)

        else:
            for iter_freq_n in range(len(freq)):
                for iter_freq_m in range(len(freq)):
                    freq_n = freq[iter_freq_n]
                    freq_m = freq[iter_freq_m]
                    out_M_im_temp[iter_freq_n, iter_freq_m] = inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type)
                    
    #-------------------------------------------------------------------------------------------------------------------           
       
    if der_used == "2nd-order":
        if std_freq/mean_freq < 1:  #%(error in frequency difference <1  #% make sure that the terms are evenly distributed)
            for iter_freq_n in range(len(freq)):
                freq_n = freq[iter_freq_n]
                freq_m = freq[0]
                C[iter_freq_n, 0] = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type)
                    
            for iter_freq_m in range(len(freq)):

                freq_n = freq[0]
                freq_m = freq[iter_freq_m]
                R[0, iter_freq_m] = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type)
                
            out_M_im_temp = toeplitz(C,R)

        else: #%if log of tau is not evenly distributed
            for iter_freq_n in range(len(freq)):
                for iter_freq_m in range(len(freq)):
                    freq_n = freq[iter_freq_n]
                    freq_m = freq[iter_freq_m]
                    out_M_im_temp[iter_freq_n, iter_freq_m] = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type)

    out_M_im[2::, 2::] = out_M_im_temp
    
    return out_M_im

##############################################################################
#----------------------------------------------------------------------------#
##############################################################################
    
def compute_A_im(freq,L):
    omega  = 2*np.pi*freq
    tau    = 1.0/freq
    N_freq = len(freq)
    
    out_A_im = np.zeros((N_freq,N_freq+2))
    if L == 1:
        out_A_im[::,0] = -2*np.pi*freq[::]
    
    for p in range(N_freq):
        for q in range(N_freq):
            if  q ==0:
                out_A_im[p, q+2] = 0.5*((omega[p]*tau[q])/(1+(omega[p]*tau[q])**2)) *np.log(tau[q+1]/tau[q])
            elif q == N_freq-1:
                out_A_im[p, q+2] = 0.5*((omega[p]*tau[q])/(1+(omega[p]*tau[q])**2))*np.log(tau[q]/tau[q-1])
            else:
                out_A_im[p, q+2] = 0.5*((omega[p]*tau[q])/(1+(omega[p]*tau[q])**2))*np.log(tau[q+1]/tau[q-1])
    
    return out_A_im

def compute_A_re(freq):
    omega  = 2*np.pi*freq
    tau    = 1.0/freq
    N_freq = len(freq)
    
    out_A_re = np.zeros((N_freq,N_freq+2))
    out_A_re[::,1] = 1
    
    for p in range(N_freq):
        for q in range(N_freq):
            if  q ==0:
                out_A_re[p, q+2] = 0.5*((1.0)/(1+(omega[p]*tau[q])**2)) *np.log(tau[q+1]/tau[q])
            elif q == N_freq-1:
                out_A_re[p, q+2] = 0.5*((1.0)/(1+(omega[p]*tau[q])**2))*np.log(tau[q]/tau[q-1])
            else:
                out_A_re[p, q+2] = 0.5*((1.0)/(1+(omega[p]*tau[q])**2))*np.log(tau[q+1]/tau[q-1])
    return out_A_re

#=====================================================================================================================
#translation of all DRT functions from matlab finished. The main-window will be defined as separate function
#=====================================================================================================================