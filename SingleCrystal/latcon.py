# latcon.py
# Program to refine just the unit cell lattice constants but not the 
# orientation matrix.
from __future__ import print_function

# A. Schultz
# January 2015

# Crystal system:
# 1 = Triclinic
# 2 = Monoclinic, b-unique
# 3 = Orthorhombic
# 4 = Tetragonal
# 5 = Rhombohedral
# 6 = Hexagonal
# 7 = Cubic

if sys.version_info > (3,):
    from tkinter import *
    import tkinter.simpledialog as tkSimpleDialog
else:
    from Tkinter import *
    import tkSimpleDialog

import sys
import numpy
from scipy.optimize import leastsq

class MyDialog(tkSimpleDialog.Dialog):

    def body(self, master):
    
        Message(master, 
            text = (
                "Initial values are obtained from latcon.inp if it exists."
                + "\nNew latcon.inp and latcon.out files will be written."
                + "\nCrystal class:"
                + "\n1 = Triclinic"
                + "\n2 = Monoclinic, b-unique"
                + "\n3 = Orthorhombic"
                + "\n4 = Tetragonal"
                + "\n5 = Rhombohedral"
                + "\n6 = Hexagonal"
                + "\n7 = Cubic"),
            aspect = 200, background = 'yellow').grid(row=0, columnspan=2)
            
        # Read user input starting parameters
        parameters = []
        try:
            user_input = open('latcon.inp', 'r')
            while True:
                lineString = user_input.readline()
                lineList = lineString.split()
                if len( lineList ) == 0: break
                if lineList[0] == '#': continue
                parameters.append(lineList[0])
            user_input.flush()
            user_input.close()
        except IOError:
            for i in range(11):
                parameters.append('')
                
        self.label_text = []
        self.label_text.append( "Crystal class: " )
        self.label_text.append( "a: " )
        self.label_text.append( "b: " )
        self.label_text.append( "c: " )
        self.label_text.append( "alpha: " )
        self.label_text.append( "beta: " )
        self.label_text.append( "gamma: " )
        self.label_text.append( "ISAW peaks or integrate file name: " )
        self.label_text.append( "Minimum I/sigI (0 for no test): " )
        self.label_text.append( "Minimum IPK (0 for no test): " )
        self.label_text.append( "Minimum d-spacing: " )

        for i in range( len( self.label_text) ):
            j = i+1
            Label(master, text = self.label_text[i]).grid(row = j, column = 0, sticky = E)
        
        # Populate entry boxes with default values
        self.entry = []
        number_of_params = len( parameters )
        for i in range(number_of_params):
            self.entry.append( Entry( master, width = 30) )
            self.entry[i].insert( 0, parameters[i] )
            j = i+1
            self.entry[i].grid( row = j, column = 1, sticky = W )
               
        
    def apply(self):
    
        self.result = []
        for i in range ( len(self.entry) ):
            self.result.append( self.entry[i].get() )

            
def get_errors(popt, pcov):
    """ Calculate esd's from covariance matrix."""
    num_variables = len(popt)

    s_sq = (infodict['fvec']**2).sum()/(len(dsp_obs)-len(popt))  # Reduced chi squared
    pcov = pcov * s_sq
    error = [] 
    for i in range(len(popt)):
        try:
          error.append( numpy.absolute(pcov[i][i])**0.5 )
        except:
          error.append( 0.00 )
    return error

#------------------ Begin ------------------

root = Tk()
root.withdraw()
root.title("latcon input")
d = MyDialog(root) 

output = open( 'latcon.out', 'w' )
output.write( 'User input:\n' )

# Get user input
crystal_system = d.result[0]
a0 = float( d.result[1] )
b0 = float( d.result[2] )
c0 = float( d.result[3] )
alpha0 = float( d.result[4] )
beta0 = float( d.result[5] )
gamma0 = float( d.result[6] )
filename = d.result[7]
min_I_sigI = float( d.result[8] )
ipkmin = int( d.result[9] )
dmin = float( d.result[10] )

# Write or over-write latcon.inp file
user_input = open( 'latcon.inp', 'w' )
user_input.write( 
        "# Crystal class:"
        + "\n# 1 = Triclinic"
        + "\n# 2 = Monoclinic, b-unique"
        + "\n# 3 = Orthorhombic"
        + "\n# 4 = Tetragonal"
        + "\n# 5 = Rhombohedral"
        + "\n# 6 = Hexagonal"
        + "\n# 7 = Cubic"
        + "\n#\n"        )
output.write( 
        "# Crystal class:"
        + "\n# 1 = Triclinic"
        + "\n# 2 = Monoclinic, b-unique"
        + "\n# 3 = Orthorhombic"
        + "\n# 4 = Tetragonal"
        + "\n# 5 = Rhombohedral"
        + "\n# 6 = Hexagonal"
        + "\n# 7 = Cubic"
        + "\n#\n"        )

for i in range( len( d.result ) ):
    user_input.write( d.result[i] + '     # ' + d.label_text[i] + '\n' )
    output.write( d.result[i] + '     # ' + d.label_text[i] + '\n' )


output.write('\n----------------------------')
output.write('\nResults:\n')

h_array = []
k_array = []
l_array = []
dsp_obs = []
# Begin reading the peaks.
input = open( filename, 'r' )
while True:

    lineString = input.readline()
    lineList = lineString.split()
    if len(lineList) == 0: break
            
    if lineList[0] != '3':
        continue
    
    seqnum = int(lineList[1])
    h = float(lineList[2])
    k = float(lineList[3])
    l = float(lineList[4])
    col = float(lineList[5])
    row = float(lineList[6])
    chan = float(lineList[7])
    L2 = float(lineList[8])
    two_theta = float(lineList[9])
    az = float(lineList[10])
    wl = float(lineList[11])
    dsp = float(lineList[12])
    ipk = int(lineList[13])
    intI = float(lineList[14])
    sigI = float(lineList[15])
    rflg = int(lineList[16])
    
    if h == k == l == 0: continue
    if dsp < dmin: continue
    
    if min_I_sigI > 0.0:
        if intI == 0.00: continue
        if sigI == 0.00: continue
        if intI/sigI < min_I_sigI: continue
        
    if ipkmin > 0:
        if ipk < ipkmin: continue
    
    h_array.append(h)
    l_array.append(l)
    k_array.append(k)
    dsp_obs.append(dsp)

print('\nNumber of peaks =', len(h_array))
output.write( '\nNumber of peaks = %d\n' % len(h_array) )
    
h_array = numpy.array(h_array)
k_array = numpy.array(k_array)
l_array = numpy.array(l_array)
dsp_obs = numpy.array(dsp_obs)

#-----------------------------------------------        
if crystal_system == '1':           # triclinic
    output.write( '\nThe crystal system is triclinic.\n' )
    
    def d_calc_triclinic(h, k, l, a, b, c, alpha, beta, gamma):
    
        alpha_rad = numpy.radians(alpha)
        sin_alpha = numpy.sin(alpha_rad)
        cos_alpha = numpy.cos(alpha_rad)
        
        beta_rad = numpy.radians(beta)
        sin_beta = numpy.sin(beta_rad)
        cos_beta = numpy.cos(beta_rad)
        
        gamma_rad = numpy.radians(gamma)
        sin_gamma = numpy.sin(gamma_rad)
        cos_gamma = numpy.cos(gamma_rad)
        
        S11 = b**2 * c**2 * sin_alpha**2
        S22 = a**2 * c**2 * sin_beta**2
        S33 = a**2 * b**2 * sin_gamma**2
        S12 = a * b * c**2 * (cos_alpha * cos_beta - cos_gamma)
        S23 = a**2 * b * c * (cos_beta * cos_gamma - cos_alpha)
        S13 = a * b**2 * c * (cos_gamma * cos_alpha - cos_beta)
        
        V = volume_triclinic(a, b, c, alpha, beta, gamma)
    
        one_over_dsq = (1/V**2) * (S11*h**2 + S22*k**2 + S33*l**2 
            + 2.0*S12*h*k + 2.0*S23*k*l + 2.0*S13*h*l)
        d_calc = numpy.sqrt(1.0 / one_over_dsq)
        return d_calc
        
    def residuals_triclinic(p0, d_obs, h, k, l):
        a, b, c, alpha, beta, gamma = p0
        d_calc = d_calc_triclinic(h, k, l, a, b, c, alpha, beta, gamma)
        err = d_obs - d_calc
        return err
        
    def volume_triclinic(a, b, c, alpha, beta, gamma):
        alpha_rad = numpy.radians(alpha)
        sin_alpha = numpy.sin(alpha_rad)
        cos_alpha = numpy.cos(alpha_rad)
        
        beta_rad = numpy.radians(beta)
        sin_beta = numpy.sin(beta_rad)
        cos_beta = numpy.cos(beta_rad)
        
        gamma_rad = numpy.radians(gamma)
        sin_gamma = numpy.sin(gamma_rad)
        cos_gamma = numpy.cos(gamma_rad)
        
        V = a*b*c*(1.0 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 
            + 2.0*cos_alpha*cos_beta*cos_gamma)**0.5
            
        return V
        
    # popt is an array of the optimized parameters
    # pcov is the covariance matrix
    p0 = numpy.zeros(6)                  # initial values of parameters
    p0[0] = a0                           
    p0[1] = b0                           
    p0[2] = c0
    p0[3] = alpha0
    p0[4] = beta0
    p0[5] = gamma0
    
    popt, pcov, infodict, errmsg, succes = leastsq(residuals_triclinic, p0, 
        args=(dsp_obs, h_array, k_array, l_array), full_output = 1)
        
    error = get_errors(popt, pcov)  # get esd's for parameters

    a = popt[0]
    b = popt[1]
    c = popt[2]
    alpha = popt[3]
    beta = popt[4]
    gamma = popt[5]
    V = volume_triclinic(a, b, c, alpha, beta, gamma)
    
    print('\nLattice constants:')
    print(7*'  %10.6f' % (a,b,c,alpha,beta,gamma,V))

    alpha1 = alpha - 0.5*error[3]
    Va1 = volume_triclinic(a, b, c, alpha1, beta, gamma)
    alpha2 = alpha + 0.5*error[3]
    Va2 = volume_triclinic(a, b, c, alpha2, beta, gamma)
    delta_V_alpha = Va2 - Va1
    
    beta1 = beta - 0.5*error[4]
    Va1 = volume_triclinic(a, b, c, alpha, beta1, gamma)
    beta2 = beta + 0.5*error[4]
    Va2 = volume_triclinic(a, b, c, alpha, beta2, gamma)
    delta_V_beta = Va2 - Va1
    
    gamma1 = gamma - 0.5*error[5]
    Va1 = volume_triclinic(a, b, c, alpha, beta, gamma1)
    gamma2 = gamma + 0.5*error[5]
    Va2 = volume_triclinic(a, b, c, alpha, beta, gamma2)
    delta_V_gamma = Va2 - Va1
    
    sigV = V * ( (error[0]/a)**2 + (error[1]/b)**2 + (error[2]/c)**2
            + (delta_V_alpha/V)**2
            + (delta_V_beta/V)**2
            + (delta_V_gamma/V)**2            
            )**0.5
    print(7*'  %10.6f' % (error[0], error[1], error[2], error[3], 
        error[4], error[5], sigV))
        
    print('\nLattice constants:')
    print(7*'  %10.6f' % (a,b,c,alpha,beta,gamma,V))
    print(7*'  %10.6f' % (error[0], error[1], error[2], error[3], error[4], error[5], sigV))
    output.write( '\nLattice constants:\n' )
    output.write( 7*'  %10.6f' % (a,b,c,alpha,beta,gamma,V) )
    output.write( '\n' )
    output.write( 7*'  %10.6f' % (error[0], error[1], error[2], error[3], error[4], error[5], sigV) )
    
#-----------------------------------------------        
if crystal_system == '2':           # monoclinic
    output.write( '\nThe crystal system is monoclinic.\n' )

    def d_calc_monoclinic(h, k, l, a, b, c, beta):
        beta_rad = numpy.radians(beta)
        sin_sq_beta = numpy.sin(beta_rad)**2
        
        term1 = (h**2)/( (a**2) * sin_sq_beta )
        term2 = (k**2)/b**2
        term3 = (l**2)/( (c**2) * sin_sq_beta)
        term4 = 2 * h * l * numpy.cos(beta_rad) / (a * c * sin_sq_beta)
        one_over_dsq = term1 + term2 + term3 - term4
        d_calc = numpy.sqrt( 1.0 / one_over_dsq )
        return d_calc
        
    def residuals_monoclinic(p0, d_obs, h, k, l):
        a, b, c, beta = p0
        d_calc = d_calc_monoclinic(h, k, l, a, b, c, beta)
        err = d_obs - d_calc
        return err
    
    # popt is an array of the optimized parameters
    # pcov is the covariance matrix
    p0 = numpy.zeros(4)                  # initial values of parameters
    p0[0] = a0                           
    p0[1] = b0                           
    p0[2] = c0
    p0[3] = beta0
    
    popt, pcov, infodict, errmsg, succes = leastsq(residuals_monoclinic, p0, 
        args=(dsp_obs, h_array, k_array, l_array), full_output = 1)

    error = get_errors(popt, pcov)  # get esd's for parameters
    
    a = popt[0]
    b = popt[1]
    c = popt[2]
    alpha = 90
    beta = popt[3]
    beta_rad = numpy.radians(beta)
    sin_beta = numpy.sin( beta_rad )
    gamma = 90
    zero = 0.0
    V = a * b * c * numpy.sin( beta_rad)
    
    sin2 = numpy.sin( numpy.radians( beta + 0.5*error[3] ) )
    sin1 = numpy.sin( numpy.radians( beta - 0.5*error[3] ) )
    delta_V = a*b*c*( numpy.sin(sin2) - numpy.sin(sin1))
    sigV = V * ( (error[0]/a)**2 + (error[1]/b)**2 + (error[2]/c)**2 
        + (delta_V / V)**2 )**0.5
    print('\nLattice constants:')
    print(7*'  %10.6f' % (a,b,c,alpha,beta,gamma,V))
    print(7*'  %10.6f' % (error[0], error[1], error[2], zero, error[3], zero, sigV))
    output.write( '\nLattice constants:\n' )
    output.write( 7*'  %10.6f' % (a,b,c,alpha,beta,gamma,V) )
    output.write( '\n' )
    output.write( 7*'  %10.6f' % (error[0], error[1], error[2], zero, error[3], zero, sigV) )
    
#-----------------------------------------------        
if crystal_system == '3':           # orthorhombic
    output.write( '\nThe crystal system is orthorhombic.\n' )
    
    def d_calc_orthorhombic(h, k, l, a, b, c):
        one_over_dsq = ( (h**2/a**2) + (k**2/b**2) + (l**2/c**2) )
        d_calc = numpy.sqrt( 1.0 / one_over_dsq )
        return d_calc
        
    def residuals_orthorhombic(p0, d_obs, h, k, l):
        a, b, c = p0
        d_calc = d_calc_orthorhombic(h, k, l, a, b, c)
        err = d_obs - d_calc
        return err

    p0 = numpy.zeros(3)                  # initial values of parameters
    p0[0] = a0                           
    p0[1] = b0                           
    p0[2] = c0
    
    popt, pcov, infodict, errmsg, succes = leastsq(residuals_orthorhombic, p0, 
        args=(dsp_obs, h_array, k_array, l_array), full_output = 1)

    error = get_errors(popt, pcov)  # get esd's for parameters
    
    a = popt[0]
    b = popt[1]
    c = popt[2]
    alpha = beta = gamma = 90
    V = a*b*c
    sigV = V * ( (error[0]/a)**2 + (error[1]/b)**2 + (error[2]/c)**2 )**0.5
    zero = 0.0
    print('\nLattice constants:')
    print(7*'  %10.6f' % (a,b,c,alpha,beta,gamma,V))
    print(7*'  %10.6f' % (error[0], error[1], error[2], zero, zero, zero, sigV))
    output.write( '\nLattice constants:\n' )
    output.write( 7*'  %10.6f' % (a,b,c,alpha,beta,gamma,V) )
    output.write( '\n' )
    output.write( 7*'  %10.6f' % (error[0], error[1], error[2], zero, zero, zero, sigV) )
    output.write( '\n' )

#-----------------------------------------------    
if crystal_system == '4':           # tetragonal
    output.write( '\nThe crystal system is tetragonal.\n' )
    
    def d_calc_tetragonal(h, k, l, a, c):
        one_over_dsq = ( (h**2 + k**2)/a**2 ) + l**2/c**2
        d_calc = numpy.sqrt( 1.0 / one_over_dsq )
        return d_calc
        
    def residuals_tetragonal(p0, d_obs, h, k, l):
        a, c = p0
        d_calc = d_calc_tetragonal(h, k, l, a, c)
        err = d_obs - d_calc
        return err
    
    # popt is an array of the optimized parameters
    # pcov is the covariance matrix
    p0 = numpy.zeros(2)                  # initial values of parameters
    p0[0] = a0                           
    p0[1] = c0                           
    
    popt, pcov, infodict, errmsg, succes = leastsq(residuals_tetragonal, p0, 
        args=(dsp_obs, h_array, k_array, l_array), full_output = 1)
        
    error = get_errors(popt, pcov)  # get esd's for parameters
        
    a = b = popt[0]
    c = popt[1]
    alpha = beta = gamma = 90
    V = a*b*c
    sigV = V * ( 2.0*(error[0]/a)**2 + (error[1]/c)**2 )**0.5
    zero = 0.0
    print('\nLattice constants:')
    print(7*'  %10.6f' % (a,b,c,alpha,beta,gamma,V))
    print(7*'  %10.6f' % (error[0], zero, error[1], zero, zero, zero, sigV))
    output.write( '\nLattice constants:\n' )
    output.write( 7*'  %10.6f' % (a,b,c,alpha,beta,gamma,V) )
    output.write( '\n' )
    output.write( 7*'  %10.6f' % (error[0], zero, error[1], zero, zero, zero, sigV) )
    
#-----------------------------------------------        
if crystal_system == '5':           # rhombohedral
    output.write( '\nThe crystal system is rhombohedral.\n' )
    
    def d_calc_rhombohedral(h, k, l, a, alpha):
        alpha_rad = numpy.radians(alpha)
        sin_alpha = numpy.radians(alpha_rad)
        cos_alpha = numpy.radians(alpha_rad)
        one_over_dsq = ( ( (h**2 + k**2 + l**2) * sin_alpha**2
            + 2.0*(h*k + k*l + h*l) * (cos_alpha**2 - cos_alpha) )
            / ( a**2 * ( 1.0 - 3.0*cos_alpha**2 + 2.0*cos_alpha**3) ) )            
        d_calc = numpy.sqrt( 1.0 / one_over_dsq )
        return d_calc
        
    def residuals_rhombohedral(p0, d_obs, h, k, l):
        a, c = p0
        d_calc = d_calc_rhombohedral(h, k, l, a, alpha)
        err = d_obs - d_calc
        return err
        
    # popt is an array of the optimized parameters
    # pcov is the covariance matrix
    p0 = numpy.zeros(2)                  # initial values of parameters
    p0[0] = a0
    p0[1] = alpha0    
    
    popt, pcov, infodict, errmsg, succes = leastsq(residuals_rhombohedral, p0, 
        args=(dsp_obs, h_array, k_array, l_array), full_output = 1)
        
    error = get_errors(popt, pcov)  # get esd's for parameters
    
    a = b = c = popt[0]
    alpha = beta = gamma = popt[1]
    cos_alpha = numpy.cos( numpy.radians(alpha) )
    V = a**3 * (1.0 - 3.0*cos_alpha**2 + 2.0*numpy.power(cos_alpha, 3.0))
    # This is not correct. Needs work.
    sigV = V * ( 3.0*(error[0]/a)**2 + 5.0*(error[1]/alpha)**2 )**0.5
    zero = 0.0
    print('\nLattice constants:')
    print(7*'  %10.6f' % (a,b,c,alpha,beta,gamma,V))
    print(7*'  %10.6f' % (error[0], zero, error[1], zero, zero, zero, sigV))
    
#-----------------------------------------------        
if crystal_system == '6':           # hexagonal
    output.write( '\nThe crystal system is hexagonal.\n' )
    
    def d_calc_hexagonal(h, k, l, a, c):
        one_over_dsq = (4.0/3.0) * ((h**2 + h*k + k**2) / a**2 ) + l**2/c**2
        d_calc = numpy.sqrt( 1.0 / one_over_dsq )
        return d_calc
        
    def residuals_hexagonal(p0, d_obs, h, k, l):
        a, c = p0
        d_calc = d_calc_hexagonal(h, k, l, a, c)
        err = d_obs - d_calc
        return err

    # popt is an array of the optimized parameters
    # pcov is the covariance matrix
    p0 = numpy.zeros(2)                  # initial values of parameters
    p0[0] = a0
    p0[1] = c0    
    
    popt, pcov, infodict, errmsg, succes = leastsq(residuals_hexagonal, p0, 
        args=(dsp_obs, h_array, k_array, l_array), full_output = 1)
        
    error = get_errors(popt, pcov)  # get esd's for parameters
    
    a = b = popt[0]
    c = popt[1]
    alpha = beta = 90
    gamma = 120
    V = numpy.sqrt(3.0) * a**2 * c / 2.0
    sigV = V * ( 2.0*(error[0]/a)**2 + (error[1]/c)**2 )**0.5
    zero = 0.0
    print('\nLattice constants:')
    print(7*'  %10.6f' % (a,b,c,alpha,beta,gamma,V))
    print(7*'  %10.6f' % (error[0], zero, error[1], zero, zero, zero, sigV))
    output.write( '\nLattice constants:\n' )
    output.write( 7*'  %10.6f' % (a,b,c,alpha,beta,gamma,V) )
    output.write( '\n' )
    output.write( 7*'  %10.6f' % (error[0], zero, error[1], zero, zero, zero, sigV) )

#-----------------------------------------------        
if crystal_system == '7':           # cubic
    output.write( '\nThe crystal system is cubic.\n' )
    
    def d_calc_cubic(h, k, l, a):
        one_over_dsq = (h**2 + k**2 + l**2) / a**2   # 1/d^2
        d_calc = numpy.sqrt( 1.0 / one_over_dsq )
        return d_calc
        
    def residuals_cubic(a, d_obs, h, k, l):
        d_calc = d_calc_cubic(h, k, l, a)
        err = d_obs - d_calc
        return err
        
    # popt is an array of the optimized parameters
    # pcov is the covariance matrix
    p0 = numpy.zeros(1)                  # initial values of parameters
    p0[0] = a0                           
    
    popt, pcov, infodict, errmsg, succes = leastsq(residuals_cubic, p0, 
        args=(dsp_obs, h_array, k_array, l_array), full_output = 1)
        
    error = get_errors(popt, pcov)  # get esd's for parameters
    
    a = b = c = popt[0]
    alpha = beta = gamma = 90
    V = a*b*c
    sigV = V * ( 3.0*(error[0]/a)**2 )**0.5
    zero = 0.0
    print('\nLattice constants:')
    print(7*'  %10.6f' % (a,b,c,alpha,beta,gamma,V))
    print(7*'  %10.6f' % (error[0], zero, zero, zero, zero, zero, sigV))
    output.write( '\nLattice constants:\n' )
    output.write( 7*'  %10.6f' % (a,b,c,alpha,beta,gamma,V) )
    output.write( '\n' )
    output.write( 7*'  %10.6f' % (error[0], zero, zero, zero, zero, zero, sigV) )
    
########################################################################    

output.write( '\n' )

print('\nResults saved to latcon.out file.')
print('\nAll done!')    

    
    

    
    

