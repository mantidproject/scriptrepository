from __future__ import print_function
import time as time
import math
import numpy
from scipy.interpolate import interp1d
from mantid.simpleapi import *
try:
  from mantidplot import *
except ImportError:
  pass
'''
chop(inst,ei,chop_type,frequency):
python implementation of CHOP ver 1.0
simulates flux and resolution for the chopper instruments
MAPS MARI and HET
original FORTRAN by TGP
matlab version JWT


Chop can either be run from the command line in matlab
chop(inst,ei,chop_type,frequency)
either
chop('mari',single or vector,'s',single or vector for frequency)
or
chop('mari',single or vector for energy,'s','12s')
requires
ei = incident energy
chop_type -- type of chopper
frequency
inst type
'''
    #first some instrument parameters
    #Instument details

global x0, xa, x1, x2, wa_mm, ha_mm, wa, ha, pslit
global dslat, radius, rho, tjit, mod_type, s, thetam, mod_type
global imod, sx, sy, sz, isam, gam, ia, ix, idet, dd, tbin, titledata


def calc_chop(ei,frequency,type):
    global x0, xa, x1, x2, wa_mm, ha_mm, wa, ha, pslit
    global dslat, radius, rho, tjit, mod_type, s, thetam, mod_type
    global imod, sx, sy, sz, isam, gam, ia, ix, idet, dd, tbin, titledata, instname, chop_par, chop_type

    #Set all of the variables required for the calculation inside this subroutine, in order to allow this single func to be called by a fitting algorithm
    #Previously, if we wished to change the instrument parameters, such as the chopper slit thickness say, we would have to re-initialise the inst, which 
    #would not work if iteratively fitting such a parameter with a single func.


########################################################
#Start of section where global instrument parameters are set
########################################################
    #For HET pseudo merlin
    x0 = 10.00
    xa = 7.190
    x1 = 1.820
    x2 = 2.5
    #NB - Helen to provide more accurate estimate of beamline apertures wa and ha
    wa_mm = 66.6670
    ha_mm = 66.6670
    wa = ha_mm / 1000.00
    ha = ha_mm / 1000.00

    # chopper details
    # now some moderator details
    # for 300K H2O
    s=numpy.zeros(6)
    #Value of s[1] from CHOP (HET) is 38.6. Use 80 to approx double the pulse width. May need to be altered after detailed comparison to data...
    s[1] = 80.00
    s[2] = 0.52260
    s[3] = 0.00
    s[4] = 0.00
    s[5] = 0.00
    th_deg = 26.70
    imod = 2
    mod_type = 'AP'

    # sample details
    sx_mm = 2.00
    sy_mm = 40.00
    sz_mm = 40.00
    isam = 0
    gam_deg = 0.00
    ia = 0
    ix = 0

    # detector details
    idet    = 2
    dd_mm   = 25.0
    tbin_us = 0.00
    
    thetam = th_deg*(math.pi/180.00)
    # function sigset set a common variable in this case to zero
    # sigset (0.0d0, 0.0d0, 0.0d0, 0.0d0)
    sx = sx_mm / 1000.00
    sy = sy_mm / 1000.00
    sz = sz_mm / 1000.00
    gam = gam_deg*math.pi/180.00

    dd = dd_mm / 1000.00
    tbin = tbin_us * 1.0e-6

    error_message = 'Chopper type is not supported for '
    
    chop_type=type
    instname='mer'

    #
    #Check whether a pre-defined chopper has been specified (using a single letter string, such as 's', 'a' (see below)
    #Or alternatively we can supply a list of chopper parameters, chop_par
    if isinstance(chop_type, (list,tuple)):
	if len(chop_type)!=6:
	    emess='Wrong number of chopper parameters specified'
	    print(emess)
	    return(emess,1)
	else:
	    chop_par=chop_type
        flux_fudge=3
    else:
		err = 0
		if type == 'c':
			chop_par=[ 1.710, 0.550, 49.00,  580.00, 0.00, 0.00 ]
			flux_fudge=3
			#print 'HET C (100meV) chopper chosen'
			titledata='HET C (100meV)'
		elif type == 'd':
			chop_par=[ 1.520, 0.550, 49.00,  410.00, 0.00, 0.00 ]
			flux_fudge=3
			#print 'HET D (50meV) chopper chosen'
			titledata='HET D (50meV)'
		elif type == 's':
			chop_par=[2.280, 0.550, 49.00, 1300.00, 0.00, 0.00 ]
			flux_fudge=3
			#print 'HET S (sloppy) chopper chosen'
			titledata='HET S (sloppy)'
		elif type == 'b':
			chop_par= [1.290, 0.550, 49.00,  920.00, 0.00, 0.00 ]
			flux_fudge=3
			#print 'HET B (200meV) chopper chosen'
			titledata='HET B (200meV)'
		elif type == 'a':
			chop_par= [0.760, 0.550, 49.00, 1300.00, 0.00, 0.0 ]
			flux_fudge=3
			#print 'HET A (500meV) chopper chosen'
			titledata='HET A (500meV)'
		elif type=='g':
			chop_par=[0.2, 0.02, 5.00, 1000000.00, 0.00, 0.00 ]
			flux_fudge=3
			#print 'MERLIN G (gadolinium) chopper chosen'
			titledata='MERLIN G'
		else:
			err=error_message+instname
			print(err)
			return (err,1)

    # Convert instrument parameters for the program (set as globals)
    pslit  = chop_par[0] / 1000.00
    dslat  = (chop_par[0] + chop_par[1]) / 1000.00
    radius = chop_par[2] / 1000.00
    rho    = chop_par[3] / 1000.00
    omega = frequency*(2*math.pi)
    tjit   = chop_par[5] * 1.0e-6

########################################################
#End of section where global instrument parameters are set
########################################################

########################################################
#Section where some other (non-instrument) parameters are set
########################################################
    fac=.95
    en_lo=0.05*ei
    eps_min=en_lo
    eps_max=fac*ei +(1-fac)*en_lo
    en_hi=eps_max

    convert = 2.3548200
    
    ###########################
    #End of extra pars section
    ###########################
    
    # calcs flux and resoltion over a range
    # en_lo some lower range for energy for the calculation
    v_lo = max( 437.391580*math.sqrt((en_lo)), 2.00*omega/(1.00/rho + (2.00*pslit/radius**2)) )
    g_lo = max( -4.00, (2.00*radius**2/pslit)*(1.00/rho - 2.00*omega/v_lo) )
    v_hi = 437.391580*math.sqrt((en_hi))
    g_hi = min( 4.00, (2.00*radius**2/pslit)*(1.00/rho - 2.00*omega/v_hi) )
    #now get the flux and the resolution width

    if numpy.size(omega)==1:
        flux=sam_flux(ei,omega)
        #Adjust by empirical flux fudge factor
        flux=flux/flux_fudge
        v_van,tsqmod,tsqchop=van_var(ei,omega)
        #do the conversions for the resolutions
        #flux = numpy.real(flux)
        van_el = numpy.real( convert * 8.747832e-4 * math.sqrt(ei**3) * ( math.sqrt(v_van + (tbin**2/12.00)) * 1.0e6 ) / x2 )

	#Put in pre-factors for tsqmod and tsqchop to allow comparison to CHOP, by giving time width at the detectors
    #Only consider elastic case, (though inelastic is calculable by setting r=(vi/vf)^3 = (Ei/Ef)^1.5
    r=1
    tsqmod_sam=(((x1+(r*x2))/x0)**2) * tsqmod
    tmod=numpy.sqrt(tsqmod_sam)
    #
    tsqchop_sam=(1+ (x1+(r*x2))/x0)**2 * tsqchop
    tchp=numpy.sqrt(tsqchop_sam)
    
    
    #eeps=range(int(eps_min),int(eps_max),1)
    eeps=numpy.linspace(eps_min,eps_max,num=19,endpoint=True)
    van=numpy.zeros(numpy.size(eeps))
    for i in range(numpy.size(eeps)):
       etrans=ei-eeps[i]
       v_van2,tsqmod2,tsqchop2=van_var(ei,omega,etrans)
       van[i] = numpy.real( convert * 8.747832e-4 * math.sqrt((ei-etrans)**3) * ( math.sqrt(v_van2 + (tbin**2/12.00)) * 1.0e6 ) / x2 )

    return van_el,van,flux,tmod,tchp,v_van
        

def van_var(*args):
    """
    calcuate vanadium widths at ei and etrans
    vanvar(ei,omega,etrans):
    """
    global x0, xa, x1, x2, wa_mm, ha_mm, wa, ha, pslit
    global dslat, radius, rho, tjit, mod_type, s, thetam, mod_type
    global imod, sx, sy, sz, isam, gam, ia, ix, idet, dd, tbin, titledata

    if len(args) == 3:
        ei =args[0]
        omega=args[1]
        eps=args[2]

    else:
        ei =args[0]
        omega=args[1]
        eps=0



    wi = 0.69468875*math.sqrt(ei)
    wf = 0.69468875*math.sqrt(ei-eps)

#
# !  get a load of widths:
# !  ---------------------
# !  moderator:
    if imod == 0:
        tsqmod=tchi(s[1]/1000.00, ei)
    elif imod == 1:
        tsqmod=tikeda(s[1], s[2], s[3], s[4], s[5], ei)
    elif imod == 2:
        tsqmod=tchi_2((s[1]/1000.00), (s[2]/1000.00), ei)

#!  chopper:
    tsqchp,ifail=tchop(omega, ei)
    ifail
    if (ifail != 0):
        tsqchp = 0.0


#!  chopper jitter:
    tsqjit = tjit**2


    v_x, v_y, v_z=sam0(6)
    atten=1.00
    dx=0.00
    dy=0.00
    dz=0.00
    v_xy=0.00

# !  detector:
    deld, sigd, sigdz, sigdd, effic=detect2(1.00, 1.00, wf, 6)
    v_dd=sigdd**2

# !  now calculate the inelastic vanadium width:
# !  -------------------------------------------
    phi=0
    v_van=van_calc(tsqmod, tsqchp, tsqjit, v_x, v_y, v_xy, v_dd, ei, eps, phi,omega)
    return v_van,tsqmod,tsqchp


#moderator functions

def tikeda(S1,S2,B1,B2,EMOD,ei):
    global x0, xa, x1, x2, wa_mm, ha_mm, wa, ha, pslit
    global dslat, radius, rho, tjit, mod_type, s, thetam, mod_type
    global imod, sx, sy, sz, isam, gam, ia, ix, idet, dd, tbin, titledata
    SIG=math.sqrt( (S1*S1) + ((S2*S2*81.8048)/ei) )
    A = 4.37392e-4 * SIG * math.sqrt(ei)
    for j in range(len(ei)):
        if (ei[j] > 130.0):
            B[j]=B2
        else:
            B[j]=B1


        R=math.math.exp(-ei/EMOD)
        tausqr[j]=(3.0/(A*A)) + (R*(2.0-R))/(B[j]*B[j])

    # variance currently in mms**2. Convert to sec**2

    return tausqr*1.0e-12

def tchi(DELTA,ei):
    global x0, xa, x1, x2, wa_mm, ha_mm, wa, ha, pslit
    global dslat, radius, rho, tjit, mod_type, s, thetam, mod_type
    global imod, sx, sy, sz, isam, gam, ia, ix, idet, dd, tbin, titledata
    VEL=437.392*math.sqrt(ei)
    tausqr=( (DELTA/1.96)/ VEL )**2
    return tausqr

def tchi_2(DELTA_0,DELTA_G,ei):
    global x0, xa, x1, x2, wa_mm, ha_mm, wa, ha, pslit
    global dslat, radius, rho, tjit, mod_type, s, thetam, mod_type
    global imod, sx, sy, sz, isam, gam, ia, ix, idet, dd, tbin, titledata
    VEL=437.392*math.sqrt(ei)
    tausqr=( ( (DELTA_0+DELTA_G*math.sqrt(ei))/1.96) / VEL )**2
    return tausqr
# end of moderator functions

def tchop(omega,ei):
#chopper
    global x0, xa, x1, x2, wa_mm, ha_mm, wa, ha, pslit
    global dslat, radius, rho, tjit, mod_type, s, thetam, mod_type
    global imod, sx, sy, sz, isam, gam, ia, ix, idet, dd, tbin, titledata
    p=pslit
    R=radius
    rho=rho
    w=omega
    ei=ei

    if (p == 0.00 and R == 0.00 and rho == 0.00):
        ierr=1
        tausqr = 0.00


# !  Calculate parameter gam:
# ! --------------------------
    veloc=437.3920*math.sqrt(ei)
    gammm=( 2.00*(R**2)/p ) * abs(1.00/rho - 2.00*w/veloc)

#!  Find regime and calculate variance:
#! -------------------------------------
    #for j in range(len(ei)):
    groot=0
    if (gammm >= 4.00):
        ierr=1
        tausqr=0.00
    else:
        ierr=0
    if (gammm <= 1.00):
        gsqr=(1.00-(gammm**2)**2 /10.00) / (1.00-(gammm**2)/6.00)
    else:
        groot=math.sqrt(gammm)
        gsqr=0.60*gammm*((groot-2.00)**2)*(groot+8.00)/(groot+4.00)

    tausqr=( (p/(2.00*R*w))**2/ 6.00) * gsqr
    return tausqr,ierr


def sam0(iout):
    global x0, xa, x1, x2, wa_mm, ha_mm, wa, ha, pslit
    global dslat, radius, rho, tjit, mod_type, s, thetam, mod_type
    global imod, sx, sy, sz, isam, gam, ia, ix, idet, dd, tbin, titledata
    ##Old code
    #varx=1
    #vary=1
    #varz=1
    ##End of old code
    #
    ##RAE code (assumes plate sample, so correct for MAPS vanadium) more sophisticated version would do different things depending on sample type
    ##but usually this contribution is small, and in any case this will be close enough for most geometries
    varx=0
    #vary=((sx)**2 + (sy)**2) #WRONG
    vary=(sy**2)/12.00
    varz=0
    return varx, vary, varz


def detect2(wd,hd,wf,iout):
	global x0, xa, x1, x2, wa_mm, ha_mm, wa, ha, pslit
	global dslat, radius, rho, tjit, mod_type, s, thetam, mod_type
	global imod, sx, sy, sz, isam, gam, ia, ix, idet, dd, tbin, titledata
	unity=1.0
	ierr=0
	unity=1.0

	if idet==1:
		# % !  Davidson type scintillator detector
		effic=1
		delta=0
		sigd =wd/math.sqrt(12.0)
		sigdz=hd/math.sqrt(12.0)
		sigdd=math.sqrt(v_dd)
		print('Li detector not supported in Pychop')
		return delta, sigd, sigdz, sigdd
	elif idet == 2:
		# % !  He cylindrical detectors binned together
		# % ! -----------------------------------------
		rad=dd/2.0
		atms=10.0
		t2rad=0.063
		effic,delta,ddsqr,v_dd,v_d=detect_he(wf,rad,atms,t2rad)

		sigd =wd/math.sqrt(12.0)
		sigdz=hd/math.sqrt(12.0)
		sigdd=math.sqrt(v_dd)
	else:
		# % !  He cylindrical detector
		# % ! -----------------------
		rad=dd/2.0
		atms=10.0
		t2rad=0.063
		effic,delta,ddsqr,v_dd,v_d=detect_he(wf,rad,atms,t2rad)
		ndet=max( anint(wd/dd),unity)
		space=2.0*rad
		v_d=v_d + (space**2) * (ndet**2 -1)/12.0
		sigd=math.sqrt(v_d)
		sigdz=hd/math.sqrt(12.0)
		sigdd=math.sqrt(v_dd)

	return delta,sigd,sigdz,sigdd,effic


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%RAE added subroutines
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def detect_he(wvec,rad,atms,t2rad):
	# % !
	# % ! Arguments
	# % ! ----------
	# % ! Entry:
	# % !   radius (m), wave-vector(A-1), no.atmospheres He (atms),
	# % !  ratio wall thickness to radius
	# % !
	# % ! Exit:
	# % !  (The above are unchanged)
	# % !  effic     efficiency (0< eff <1)
	# % !  delta     shift in effective position of detector from centre (m)
	# % !  ddsqr     mean square depth of absorption (w.r.t. centre) (m**2)
	# % !  v_dd      variance of depth (m**2)
	# % !  v_d       variance of width (m**2)
	# % !  ierr   =  0  no problems
	# % !         =  1  rad< 0, t2rad <0 or >1, atms <0, wvec too small
	# % !         =  2  error in called routines
	# % !
	# % ! Origin of data
	# % ! --------------
	# % !  CKL data :
	# % !   "At 2200 m/s xsect=5327 barns    En=25.415 meV         "
	# % !   "At 10 atms, rho_atomic=2.688e-4,  so sigma=1.4323 cm-1"
	# % !
	# % !  These data are not quite consistent, but the errors are small :
	# % !    2200 m/s = 25.299 meV
	# % !    5327 barns & 1.4323 cm-1 ==> 10atms ofideal gas at 272.9K
	# % !   but at what temperature are the tubes "10 atms" ?
	# % !
	# % !  Shall use  1.4323 cm-1 @ 3.49416 A-1 with sigma prop. 1/v
	# % !
	# % !  This corresponds to a reference energy of 25.299meV, NOT 25.415.
	# % ! This accounts for a difference of typically 1 pt in 1000 for
	# % ! energies around a few hundred meV. The approximate quadrature
	# % ! of CKL is accurate typically to 1 pt in 10000 for efficiency.
	# % ! The Chebyshev approximation used in this routine is accurate to
	# % ! 1 pt in 10**12 for all the returned moments. The routine is
	# % ! typically 5 times faster than CKL.
	# % !

	sigref=143.23
	wref=3.49416
	atmref=10.0
	const = sigref*wref/atmref

	effic=0.0
	delta=0.0
	ddsqr=0.0
	v_dd =0.0
	v_d  =0.0
	

	if rad < 0.0 or t2rad < 0.0 or t2rad > 1.0 or atms < 0.0:
		print('Error with detect_he input parameters')
		return
	else:
		reff=( rad*(1.0-t2rad) )
		var =2.0 * ( rad*(1.0-t2rad) ) * ( const*atms )
		if wvec < (var*1.0e-18):
			error('Error with size of wavevector for input pars')
		else:
			alf=var/wvec
			[effic,delta,ddsqr,v_dd,v_d]=tube_mts(alf)
			delta=delta*reff
			ddsqr=ddsqr*(reff**2)
			v_dd =v_dd *(reff**2)
			v_d  =v_d  *(reff**2)

	return effic,delta,ddsqr,v_dd,v_d
#%============

def tube_mts(alf):

	# % !
	# % !  T.G.Perring 6/4/90
	# % !
	# % !  Given ALF (radius in m.f.p.), calculates:
	# % !    EFF   efficiency
	# % !    DEL   mean depth of absorption (w.r.t centre)
	# % !    XSQR  mean of (depth**2) (w.r.t. centre)
	# % !    VX    variance of depth of absorption
	# % !    VY    variance of width
	# % !
	# % !  The routine approximates the functions by Chebyshev polynomial
	# % ! expansions over the ranges  0 =< alf =<9  and
	# % ! 10 =< alf =< (infinity). For  9 =< alf =< 10 a linear combination
	# % ! of the two approximations is taken.
	# % !
	# % !  The routine gives relative and absolute accuracy of the quantities
	# % ! to better than 10**-12 for all positive ALF.
	# % !
	# % !    IERR  returned as  0  ALF .ge. 0.0
	# % !                       1  ALF .lt. 0.0
	# % !
	g0=(32.0-3.0*(math.pi**2))/48.0
	g1=14.0/3.0-(math.pi**2)/8.0



	c_eff_f=[0.7648360390553052,
		-0.3700950778935237, 0.1582704090813516,
		-6.0170218669705407e-02, 2.0465515957968953e-02,
		-6.2690181465706840e-03, 1.7408667184745830e-03,
		-4.4101378999425122e-04, 1.0252117967127217e-04,
		-2.1988904738111659e-05, 4.3729347905629990e-06,
		-8.0998753944849788e-07, 1.4031240949230472e-07,
		-2.2815971698619819e-08, 3.4943984983382137e-09,
		-5.0562696807254781e-10, 6.9315483353094009e-11,
		-9.0261598195695569e-12, 1.1192324844699897e-12,
		-1.3204992654891612e-13, 1.4100387524251801e-14,
		-8.6430862467068437e-16,-1.1129985821867194e-16,
		-4.5505266221823604e-16, 3.8885561437496108e-16]

	c_eff_g=[ 2.033429926215546,
		-2.3123407369310212e-02, 7.0671915734894875e-03,
		-7.5970017538257162e-04, 7.4848652541832373e-05,
		4.5642679186460588e-05,-2.3097291253000307e-05,
		1.9697221715275770e-06, 2.4115259271262346e-06,
		-7.1302220919333692e-07,-2.5124427621592282e-07,
		1.3246884875139919e-07, 3.4364196805913849e-08,
		-2.2891359549026546e-08,-6.7281240212491156e-09,
		3.8292458615085678e-09, 1.6451021034313840e-09,
		-5.5868962123284405e-10,-4.2052310689211225e-10,
		4.3217612266666094e-11, 9.9547699528024225e-11,
		1.2882834243832519e-11,-1.9103066351000564e-11,
		-7.6805495297094239e-12, 1.8568853399347773e-12]

	c_del_f=[1.457564928500728,
		-0.2741263150129247,     1.4102406058428482e-02,
		1.1868136977190956e-02,-4.7000120888695418e-03,
		6.7071002620380348e-04, 1.2315212155928235e-04,
		-8.7985748380390304e-05, 1.8952644758594150e-05,
		4.4101711646149510e-07,-1.5292393205490473e-06,
		4.5050196748941396e-07,-2.9971703975339992e-08,
		-2.3573145628841274e-08, 9.6228336343706644e-09,
		-1.3038786850216866e-09,-2.9423462000188749e-10,
		1.8813720970012326e-10,-3.7682054143672871e-11,
		-1.9125961925325896e-12, 3.3516145414580478e-12,
		-9.0842416922143343e-13, 4.3951786654616853e-14,
		4.5793924208226145e-14,-1.4916540225229369e-14]

	c_del_g=[1.980495234559052,
		1.3148750635418816e-02,-3.5137830163154959e-03,
		1.4111112411286597e-04,-2.4707009281715875e-05,
		-4.9602024972950076e-08, 1.5268651833078018e-06,
		-4.8070752083129165e-07,-3.5826648758785495e-08,
		6.0264253483044428e-08,-4.2948016776289677e-09,
		-7.5840171520624722e-09, 1.0468151234732659e-09,
		1.1267346944343615e-09,-1.4810551229871294e-10,
		-1.9605287726598419e-10, 9.8596597553068932e-12,
		3.6752354493074790e-11, 3.2634850377633029e-12,
		-6.6207839211074316e-12,-1.9158341579839089e-12,
		9.6091495871419851e-13, 6.3198529742791721e-13,
		-6.4681177081027385e-14,-1.8198241524824965e-13]

	c_xsqr_f=[2.675986138240137,
		0.4041429091631520,     2.1888771714164858e-02,
		-3.4310286472213617e-02, 9.8724790919419380e-03,
		-7.7267251256297631e-04,-4.6681418487147020e-04,
		2.0604262514245964e-04,-3.1387761886573218e-05,
		-5.1728966665387510e-06, 3.9417564710109155e-06,
		-8.6522505504893487e-07,-1.6220695979729527e-08,
		6.8546255754808882e-08,-2.0405647520593817e-08,
		1.4047699248287415e-09, 1.0523175986154598e-09,
		-4.3422357653977173e-10, 5.9649738481937220e-11,
		1.3017424915773290e-11,-8.4605289440986553e-12,
		1.7046483669069801e-12, 8.2185647176657995e-14,
		-1.4448442442471787e-13, 3.5720454372167865e-14]

	c_xsqr_g=[1.723549588238691,
		0.1365565801015080,     2.0457962179522337e-03,
		-3.9875695195008110e-04, 2.3949621855833269e-05,
		-1.6129278268772751e-06,-1.1466609509480641e-06,
		4.3086322193297555e-07, 1.7612995328875059e-09,
		-4.5839686845239313e-08, 5.9957170539526316e-09,
		5.3204258865235943e-09,-1.1050097059595032e-09,
		-7.7028480982566094e-10, 1.5644044393248180e-10,
		1.3525529252156332e-10,-1.5409274967126407e-11,
		-2.6052305868162762e-11,-8.3781981352615275e-13,
		4.8823761700234058e-12, 1.1086589979392158e-12,
		-7.5851658287717783e-13,-4.0599884565395428e-13,
		7.9971584909799275e-14, 1.3500020545897939e-13]

	c_vx_f=[1.226904583058190,
		-0.3621914072547197,     6.0117947617747081e-02,
		1.8037337764424607e-02,-1.4439005957980123e-02,
		3.8147446724517908e-03, 1.3679160269450818e-05,
		-3.7851338401354573e-04, 1.3568342238781006e-04,
		-1.3336183765173537e-05,-7.5468390663036011e-06,
		3.7919580869305580e-06,-6.4560788919254541e-07,
		-1.0509789897250599e-07, 9.0282233408123247e-08,
		-2.1598200223849062e-08,-2.6200750125049410e-10,
		1.8693270043002030e-09,-6.0097600840247623e-10,
		4.7263196689684150e-11, 3.3052446335446462e-11,
		-1.4738090470256537e-11, 2.1945176231774610e-12,
		4.7409048908875206e-13,-3.3502478569147342e-13]

	c_vx_g=[1.862646413811875,
		7.5988886169808666e-02,-8.3110620384910993e-03,
		1.1236935254690805e-03,-1.0549380723194779e-04,
		-3.8256672783453238e-05, 2.2883355513325654e-05,
		-2.4595515448511130e-06,-2.2063956882489855e-06,
		7.2331970290773207e-07, 2.2080170614557915e-07,
		-1.2957057474505262e-07,-2.9737380539129887e-08,
		2.2171316129693253e-08, 5.9127004825576534e-09,
		-3.7179338302495424e-09,-1.4794271269158443e-09,
		5.5412448241032308e-10, 3.8726354734119894e-10,
		-4.6562413924533530e-11,-9.2734525614091013e-11,
		-1.1246343578630302e-11, 1.6909724176450425e-11,
		5.6146245985821963e-12,-2.7408274955176282e-12]

	c_vy_f=[2.408884004758557,
		0.1441097208627301,    -5.0093583831079742e-02,
		1.0574012517851051e-02,-4.7245491418700381e-04,
		-5.6874753986616233e-04, 2.2050994176359695e-04,
		-3.0071128379836054e-05,-6.5175276460682774e-06,
		4.2908624511150961e-06,-8.8327783029362728e-07,
		-3.5778896608773536e-08, 7.6164115048182878e-08,
		-2.1399959173606931e-08, 1.1599700144859781e-09,
		1.2029935880786269e-09,-4.6385151497574384e-10,
		5.7945164222417134e-11, 1.5725836188806852e-11,
		-9.1953450409576476e-12, 1.7449824918358559e-12,
		1.2301937246661510e-13,-1.6739387653785798e-13,
		4.5505543777579760e-14,-4.3223757906218907e-15]

	c_vy_g=[1.970558139796674,
		1.9874189524780751e-02,-5.3520719319403742e-03,
		2.3885486654173116e-04,-4.1428357951582839e-05,
		-6.3229035418110869e-07, 2.8594609307941443e-06,
		-8.5378305322625359e-07,-8.2383358224191738e-08,
		1.1218202137786015e-07,-6.0736651874560010e-09,
		-1.4453200922748266e-08, 1.7154640064021009e-09,
		2.1673530992138979e-09,-2.4074988114186624e-10,
		-3.7678839381882767e-10, 1.1723938486696284e-11,
		7.0125182882740944e-11, 7.5127332133106960e-12,
		-1.2478237332302910e-11,-3.8880659802842388e-12,
		1.7635456983633446e-12, 1.2439449470491581e-12,
		-9.4195068411906391e-14,-3.4105815394092076e-13]


	if alf < 0:
		print('alf < 0, invalid choice')
	else:
		if alf <= 9.0:
			eff =(math.pi/4.00)*alf* chbmts(0.00,10.00,c_eff_f,25,alf)
			delta =-0.125*alf* chbmts(0.00,10.00,c_del_f,25,alf)
			xsqr=0.25* chbmts(0.00,10.00,c_xsqr_f,25,alf)
			vx  =0.25* chbmts(0.00,10.00,c_vx_f,25,alf)
			vy  =0.25* chbmts(0.00,10.00,c_vy_f,25,alf)
		elif alf >= 10.00:
			y=1.0-18.0/alf
			eff =1.00 - chbmts(-1.00,1.00,c_eff_g,25,y)/alf**2
			delta =(2.0*chbmts(-1.00,1.00,c_del_g,25,y)/alf - 0.25*math.pi) / eff
			xsqr=((-math.pi/alf)* chbmts(-1.00,1.00,c_xsqr_g,25,y) + 2.0/3.0 ) / eff
			vx  =g0 + g1*chbmts(-1.00,1.00,c_vx_g,25,y)/(alf**2)
			vy  =( -chbmts(-1.00,1.00,c_vy_g,25,y)/(alf**2) + 1.0/3.0 ) / eff
		else:
			eff_f =(math.pi/4.00)*alf*chbmts(0.00,10.00,c_eff_f,25,alf)
			del_f =-0.125*alf*chbmts(0.00,10.00,c_del_f,25,alf)
			xsqr_f=0.25*chbmts(0.00,10.00,c_xsqr_f,25,alf)
			vx_f  =0.25*chbmts(0.00,10.00,c_vx_f,25,alf)
			vy_f  =0.25*chbmts(0.00,10.00,c_vy_f,25,alf)
			y=1.0-18.0/alf
			eff_g =1.00 - chbmts(-1.00,1.00,c_eff_g,25,y)/alf**2
			del_g =(2.0*chbmts(-1.00,1.00,c_del_g,25,y)/alf - 0.25*math.pi) / eff_g
			xsqr_g=((-math.pi/alf)*chbmts(-1.00,1.00,c_xsqr_g,25,y) + 2.0/3.0 ) / eff_g
			vx_g  =g0 + g1*chbmts(-1.00,1.00,c_vx_g,25,y)/(alf**2)
			vy_g  =( -chbmts(-1.00,1.00,c_vy_g,25,y)/(alf**2) + 1.0/3.0 ) / eff_g
			eff =(10.0-alf)*eff_f  + (alf-9.0)*eff_g
			delta =(10.0-alf)*del_f  + (alf-9.0)*del_g
			xsqr=(10.0-alf)*xsqr_f + (alf-9.0)*xsqr_g
			vx  =(10.0-alf)*vx_f   + (alf-9.0)*vx_g
			vy  =(10.0-alf)*vy_f   + (alf-9.0)*vy_g


	
	return eff,delta,xsqr,vx,vy		

#%==============

def chbmts(a,b,c,m,x):

	# % !
	# % ! Essentially CHEBEV of "Numerical Recipes"
	# % !


	d=0.0
	ddd=0.0
	y=(2.0*x-a-b)/(b-a)
	y2=2.0*y

	myrange=list(range(m-1,0,-1))
	
	#print('myrange')
	#print(myrange)
	#print(c)
	
	for j in myrange:
		sv=d
		d=y2*d-ddd+c[j]
		ddd=sv

	out=y*d-ddd+0.5*c[0]

	return out


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%End of RAE added subroutines
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


def van_calc(v_mod, v_ch, v_jit, v_x, v_y, v_xy, v_dd,ei, eps, phi,omega):
#[v_van,v_van_m,v_van_ch,v_van_jit,v_van_ya,v_van_x,v_van_y,v_van_xy,v_van_dd]
    global x0, xa, x1, x2, wa_mm, ha_mm, wa, ha, pslit
    global dslat, radius, rho, tjit, mod_type, s, thetam, mod_type
    global imod, sx, sy, sz, isam, gam, ia, ix, idet, dd, tbin, titledata
    veli=437.39160*math.sqrt( ei )
    velf=437.39160*math.sqrt( ei-eps )
    rat=(veli/velf)**3


    tanthm=math.tan(thetam)

    am   = -(x1+rat*x2)/x0
    ach  = (1.00 + (x1+rat*x2)/x0)
    g1 = (1.00 - (omega*(x0+x1)*tanthm/veli))#wrong (in original CHOP!) - should be (xa+x1), not (x0+x1)
    #g1 = (1.00 - (omega*(xa+x1)*tanthm/veli))
    g2 = (1.00 - (omega*(x0-xa)*tanthm/veli) )
    f1 =  1.00 + ((x1/x0)*g1)
    f2 =  1.00 + ((x1/x0)*g2)
    gg1 = g1 / ( omega*(xa+x1) )
    gg2 = g2 / ( omega*(xa+x1) )
    ff1 = f1 / ( omega*(xa+x1) )
    ff2 = f2 / ( omega*(xa+x1) )
    aa = ( (math.cos(gam)/veli) - (math.cos(gam-phi)/velf) ) - (ff2*math.sin(gam))
    bb = ((-math.sin(gam)/veli) + (math.sin(gam-phi)/velf) ) - (ff2*math.cos(gam))
    aya  = ff1 + ((rat*x2/x0)*gg1)
    ax   = aa  - ((rat*x2/x0)*gg2*math.sin(gam))
    ay   = bb  - ((rat*x2/x0)*gg2*math.cos(gam))
    a_dd  = 1.00/velf

    v_van_m  = am**2  * v_mod
    v_van_ch = ach**2 * v_ch
    v_van_jit= ach**2 * v_jit
    v_van_ya = aya**2 * (wa**2/12.00)
    v_van_x  = ax**2  * v_x
    v_van_y  = ay**2  * v_y
    v_van_xy = ax*ay  * v_xy
    v_van_dd = a_dd**2* v_dd

    ##Old version:
    #v_van = (v_van_m + v_van_ch + v_van_jit + v_van_ya)
    
    ##New (RAE) version:
    v_van = (v_van_m + v_van_ch + v_van_jit + v_van_ya + v_van_y + v_van_dd)
    return v_van



def sam_flux(ei,omega):
    '''# Calculates the flux at the sample position for HET or MARI in  n / cm**2 . uA.s for the U target
    # (all distances in m, angles in rad, omega rad/s, ei in meV)
    '''
    global x0, xa, x1, x2, wa_mm, ha_mm, wa, ha, pslit
    global dslat, radius, rho, tjit, mod_type, s, thetam, mod_type
    global imod, sx, sy, sz, isam, gam, ia, ix, idet, dd, tbin, titledata
    flux=[]
    #flux1=flux_calc(ei)
    #print(flux1)
    #area=achop(ei,omega)
    ##for j in range(len(ei)):
    #flux_ref = 84403.060*ei*(flux1/math.cos(thetam))*(area/dslat)*(wa*ha)/(x0*(x1+xa)**2)
    #print(flux_ref)
    ##Matlab code:
    ##Measured WB flux from RIB's calibrated monitor
    #fid=fopen('C:\Russell\Instrumentation\MAPS\ResolutionCode\calibrated_flux_MAPS.dat');
    #flux_maps=fscanf(fid,'%f',[2,inf])';
    #fclose(fid);
    #flux_maps(:,1)=81.81./(flux_maps(:,1)).^2;
    
    #Python equivalent
    #maps_flux=numpy.loadtxt('/data/Instrumentation/MAPS/FluxBenchmark/calibrated_flux_MAPS.dat')
    
    #Store the info inside this file, for better portability
    maps_lam,maps_flux=flux_store()
    #Interpolate for wavelength
    f=interp1d(maps_lam,maps_flux,kind='cubic')
    
    area=achop(ei,omega)
    lam=numpy.sqrt(81.81/ei)
    flux_meas=f(lam)
    #print(flux_meas)
    flux = (84403.060)*flux_meas*(area/dslat)*1/(x0*(x1+x0))
    #print(flux)
    return flux

def flux_calc(ei):
    global x0, xa, x1, x2, wa_mm, ha_mm, wa, ha, pslit
    global dslat, radius, rho, tjit, mod_type, s, thetam, mod_type
    global imod, sx, sy, sz, isam, gam, ia, ix, idet, dd, tbin, titledata
    conv1=3.615
    conv2=9.104157e-12
    conv=conv1*conv2
    en_ev=ei/1000.00
    ch_mod=mod_type
    phi0=flux_norm(ch_mod)
    phifun=flux_fun( en_ev, ch_mod)

    flux=(conv*( phi0*math.cos(thetam) )*(math.sqrt(en_ev)*phifun))/1
    return flux

def flux_norm(ch_mod):
    global x0, xa, x1, x2, wa_mm, ha_mm, wa, ha, pslit
    global dslat, radius, rho, tjit, mod_type, s, thetam, mod_type
    global imod, sx, sy, sz, isam, gam, ia, ix, idet, dd, tbin, titledata
    if ch_mod =='A':
        phi0=1.00
    if ch_mod == 'AP':
        phi0=2.80
    if ch_mod =='H2':
        phi0=1.80
    if ch_mod == 'CH4':
        phi0=2.60

    return phi0

def flux_fun( en_ev, ch_mod):
    '''
    #%!
    #%!  Calculates the energy dependence of the flux (see Perring et al
    #%! RAL-85-029)
    #%!
    #%!  Entry
    #%!     en_ev  : eV
    #%!     ch_mod : 'A'  'AP'  'CH4'  'H2' 'TEST'
    #%!
    #%!  Exit
    #%!   phifun   : the functional dependace with energy
    '''
    global x0, xa, x1, x2, wa_mm, ha_mm, wa, ha, pslit
    global dslat, radius, rho, tjit, mod_type, s, thetam, mod_type
    global imod, sx, sy, sz, isam, gam, ia, ix, idet, dd, tbin, titledata
    if ch_mod =='A':
        ijoin=0
        rj=0.00
        t=0.00
        a=0.00
        w1=0.00
        w2=0.00
        w3=0.00
        w4=0.00
        w5=0.00
        ierr=1
    if ch_mod =='AP':
        ijoin=0
        rj=2.250
        t=0.0320
        a=0.950
        w1=120.00
        w2=10.00
        w3=0.00
        w4=0.00
        w5=0.00
    if ch_mod =='H2':
        ijoin=1
        rj=2.350
        t=0.00210
        a=0.950
        w1=15.50
        w2=3.10
        w3=11.00
        w4=0.2540
        w5=0.02750
    if ch_mod =='CH4':
        ijoin=0
        rj=2.10
        t=0.0110
        a=0.920
        w1=55.0
        w2=7.00
        w3=0.00
        w4=0.00
        w5=0.00

#
#
# ! Calculation:
# ! ------------
    phi_max=rj*(en_ev/(t**2))*math.exp(-en_ev/t)
    phi_epi=1.00/(en_ev)**a

    math.expon=math.exp(-w1/math.sqrt(1000.00*en_ev)+w2)
    delt1=math.expon/ (  1.00 + math.expon  )

    if ijoin == 1:
        math.expon=math.exp( (w4 - 1.00/math.sqrt(1000.00*en_ev) )/w5 )
        delt2=1.00 + w3/( 1.00 + math.expon )
    else:
        delt2=1.00


    phifun=phi_max + delt1*delt2*phi_epi
    return phifun
def achop(ei,omega):
    '''
    # !
    # !   Calculates the integral of the chopper transmission function
    # !  P(h,t) over time and distance for any energy. The answer is in m.S
    # !   New version 15/1/90
    # !
    # !    p       slit thickness (m)                              R*8
    # !    R       slit package diameter (m)                       R*8
    # !    rho     slit radius of curvature (m)                    R*8
    # !    w       angular frequency of rotor (rad/sec)            R*8
    # !    ei      energy the rotor has been phased for (meV)      R*8
    # !    area    intergral                                       R*8
    # !    ierr    error indicator                                 integer
    # !               =0  no problems
    # !               =1  if no transmission  AREA set to zero
    '''
    global x0, xa, x1, x2, wa_mm, ha_mm, wa, ha, pslit
    global dslat, radius, rho, tjit, mod_type, s, thetam, mod_type
    global imod, sx, sy, sz, isam, gam, ia, ix, idet, dd, tbin, titledata
    dd=dslat
    p1=pslit
    R1=radius
    rho1=rho
    w1=omega
    ei=ei

    vela=437.3920*math.sqrt(ei)
    gamm=( 2.00*(R1**2)/p1 ) * abs(1.00/rho1 - 2.00*w1/vela)

# !  Find regime and calculate variance:
# ! -------------------------------------

    #for j in range(numpy.size(ei)):
    groot=0
    if (gamm >= 4.00):
        f1=0.0
        print('no transmission at ', ei, 'meV at ',omega/(2*math.pi), 'Hz')
    else:
        ierr=0
        if gamm <= 1.00:
            f1=1.00-(gamm**2)/6.00
        else:
            groot=math.sqrt(gamm)
            f1=groot*((groot-2.00)**2)*(groot+4.00)/6.00



    area=( (p1**2)/(2.00*R1*w1) ) * f1
    return area

def frange(limit1, limit2 = None, increment = 1.):
  """
  Range function that accepts floats (and integers).

  Usage:
  frange(-2, 2, 0.1)
  frange(10)
  frange(10, increment = 0.5)

  The returned value is an iterator.  Use list(frange) for a list.
  """

  if limit2 is None:
    limit2, limit1 = limit1, 0.
  else:
    limit1 = float(limit1)

  count = int(math.ceil(limit2 - limit1)/increment)
  return (limit1 + n*increment for n in range(count))



def flux_store():
	#Return stored flux from MERLIN calibrated flux measurement. Avoids needing to keep in file with PyChop distribution

	ei=[1.5000e-01,
		2.5000e-01,
		3.5000e-01,
		4.5000e-01,
		5.5000e-01,
		6.5000e-01,
		7.5000e-01,
		8.5000e-01,
		9.5000e-01,
		1.0500e+00,
		1.1500e+00,
		1.2500e+00,
		1.3500e+00,
		1.4500e+00,
		1.5500e+00,
		1.6500e+00,
		1.7500e+00,
		1.8500e+00,
		1.9500e+00,
		2.0500e+00,
		2.1500e+00,
		2.2500e+00,
		2.3500e+00,
		2.4500e+00,
		2.5500e+00,
		2.6500e+00,
		2.7500e+00,
		2.8500e+00,
		2.9500e+00,
		3.0500e+00,
		3.1500e+00,
		3.2500e+00,
		3.3500e+00,
		3.4500e+00,
		3.5500e+00,
		3.6500e+00,
		3.7500e+00,
		3.8500e+00,
		3.9500e+00,
		4.0500e+00,
		4.1500e+00,
		4.2500e+00,
		4.3500e+00,
		4.4500e+00,
		4.5500e+00,
		4.6500e+00,
		4.7500e+00,
		4.8500e+00,
		4.9500e+00,
		5.0500e+00,
		5.1500e+00,
		5.2500e+00,
		5.3500e+00,
		5.4500e+00,
		5.5500e+00,
		5.6500e+00,
		5.7500e+00,
		5.8500e+00,
		5.9500e+00,
		6.0500e+00,
		6.1500e+00,
		6.2500e+00,
		6.3500e+00,
		6.4500e+00,
		6.5500e+00,
		6.6500e+00,
		6.7500e+00,
		6.8500e+00,
		6.9500e+00,
		7.0500e+00,
		7.1500e+00,
		7.2500e+00,
		7.3500e+00,
		7.4500e+00,
		7.5500e+00,
		7.6500e+00,
		7.7500e+00,
		7.8500e+00,
		7.9500e+00,
		8.0500e+00,
		8.1500e+00,
		8.2500e+00,
		8.3500e+00,
		8.4500e+00,
		8.5500e+00,
		8.6500e+00,
		8.7500e+00,
		8.8500e+00,
		8.9500e+00,
		9.0500e+00,
		9.1500e+00,
		9.2500e+00,
		9.3500e+00,
		9.4500e+00,
		9.5500e+00,
		9.6500e+00,
		9.7500e+00,
		9.8500e+00,
		9.9500e+00,
		1.0050e+01]
	flux=[3.4827e+07,
		1.9455e+07,
		1.3618e+07,
		1.0922e+07,
		9.1074e+06,
		8.9949e+06,
		1.0717e+07,
		1.4113e+07,
		1.7890e+07,
		2.0977e+07,
		2.2164e+07,
		2.1097e+07,
		1.9449e+07,
		1.7522e+07,
		1.5395e+07,
		1.3518e+07,
		1.1794e+07,
		1.0664e+07,
		9.4134e+06,
		7.8688e+06,
		6.7881e+06,
		5.9428e+06,
		4.9238e+06,
		4.2219e+06,
		3.8432e+06,
		3.3606e+06,
		2.9299e+06,
		2.6110e+06,
		2.3335e+06,
		2.0563e+06,
		1.8005e+06,
		1.6002e+06,
		1.4490e+06,
		1.3368e+06,
		1.2599e+06,
		1.1239e+06,
		1.0524e+06,
		9.7266e+05,
		8.6218e+05,
		7.6681e+05,
		7.0049e+05,
		6.7822e+05,
		6.2993e+05,
		6.2571e+05,
		5.8608e+05,
		5.3686e+05,
		5.5309e+05,
		5.2245e+05,
		4.8365e+05,
		4.4384e+05,
		3.9445e+05,
		3.7199e+05,
		3.2812e+05,
		3.1162e+05,
		2.9301e+05,
		2.6704e+05,
		2.6396e+05,
		2.4635e+05,
		2.1940e+05,
		2.1146e+05,
		2.0127e+05,
		1.9412e+05,
		1.7908e+05,
		1.7902e+05,
		1.6332e+05,
		1.1746e+05,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00,
		0.0000e+00]

	return ei,flux

