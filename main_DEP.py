# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 14:05:08 2017

Main file for stability mapping

@author: e.nguyen-van
"""
# import all modules
import ReadCoefAero
import AeroForces
import numpy as np
import math
import scipy.linalg
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import ATRgeometry
import equation as e
from scipy.optimize import  minimize
#import pdb
import time
#import matplotlib.image as mpimg
import sys
import MapVarClass

# --- Option de sortie console ---
np.set_printoptions(threshold=np.nan)
# fonts for graphics
font = {'family':'sherif',
        'serif':'Times New Roman',
        'size'   : 12}

plt.rc('font', **font)
# --- 

""" Algorithm set up """
Neng =12
inop_eng = 3
g = ATRgeometry.data(0.7,Neng,inop_eng) # arg = Vtsize + options(Neng, inop_eng, vertical tail parameters)
print(g.PosiEng)
# --- Test case and steady parameters
H_base = 0000 # in m the altitude
V_base = 60
beta_base = 5/180*math.pi
gamma = 2.0/180*math.pi 
R = 0 # in meters the turn radius
phimax = 5 # in degree the max bank angle authorized
alphamax=10 # in degree to adjust if it is below 71m/s
deltaRmax=0.75*30 # in degree
ThrottleMax = 1 # max thrust level
ThrottleMin = 1e-9 # min thruttle, don't accept 0 thrust
# ---- Optim parameter ------
MaxIter=100 # only for single run
tolerance=1e-5

# --- additional parameters (default edited during execution) ---
g.set_nofin(True) # =True means : no rudder used

# --- dictionnary for type of aircraft studied. aircraft: ATR72, version : 'original', 'DEPoriginal', 'DEPnofin'
g.hangar={'aircraft':'ATR72', 'version':'DEPoriginal'}

# --- plot coeff evolution
plotcoef = False
# --- plot control histogram
HistoGouv = False

# --- Study jacobian
gojac = False
storeJac = False
goFinVariation= False
FinRatioVec = np.linspace(0.1,1.1,20)
CstSpan = False
CstA = True

# --- mapping settings ---
domap = True
MapName = "Beta_Vel"      # key words "Vel", "Beta", "Gamma", "Omega", separator: "_"

MapVmax = 65              #velocity limits in degrees
MapVmin = 48
MapVstep = 1

MapBetaMax = 21
MapBetaMin = -20
MapBetaStep = 1

MapGammaMax = 5     #gamma limits in degrees
MapGammaMin = 0
MapGammaStep = 0.5

MapOmegaMax = 0.07    #omega limits
MapOmegaMin = -0.06
MapOmegaStep = 0.01

BetaStall = 15         #beta stall limit


#Flight envelop Bounding parameters
if g.nofin==True:
    BparamName=["alpha","phi","deltaa"]
    LenSmallBparam=len(BparamName)
    for i in range(Neng):
        BparamName.append("deltax"+str(i)) #complete the name
    BparamLim={'alpha':alphamax/180*math.pi,'phi':phimax/180*math.pi,'deltaa':20/180*math.pi,'deltaxmax':ThrottleMax,'deltaxmin':ThrottleMin}
    BparamThres={'alpha':0.92,'phi':0.99,'deltaa':0.95,'deltaxmax':0.98,'deltaxmin':0.01}
    LenBparam=len(BparamName)
    DispBparam={}
    for i in range(LenBparam):
        DispBparam[BparamName[i]]=np.array([]).reshape(0,2)
    BparamXvecPosi={'alpha':0,'phi':4,'deltaa':6}
    for i in range(LenSmallBparam,LenBparam):
        BparamXvecPosi[BparamName[i]]=8+i-LenSmallBparam
else:
    BparamName=["alpha","phi","deltaa","deltaR"]
    LenSmallBparam=len(BparamName)
    for i in range(Neng):
        BparamName.append("deltax"+str(i)) #complete the name
    BparamLim={'alpha':alphamax/180*math.pi,'phi':phimax/180*math.pi,'deltaa':20/180*math.pi,'deltaR':deltaRmax/180*math.pi,'deltaxmax':ThrottleMax,'deltaxmin':ThrottleMin}
    BparamThres={'alpha':0.92,'phi':0.99,'deltaa':0.95,'deltaR':0.97,'deltaxmax':0.99,'deltaxmin':0.01}
    LenBparam=len(BparamName)
    DispBparam={}
    for i in range(LenBparam):
        DispBparam[BparamName[i]]=np.array([]).reshape(0,2)
    BparamXvecPosi={'alpha':0,'phi':4,'deltaa':6,'deltaR':8}
    for i in range(LenSmallBparam,LenBparam):
        BparamXvecPosi[BparamName[i]]=9+i-LenSmallBparam
    


# --- hard coded velocity, corresponding air density and Sound vel ---
# The velocities corresponds to key points of ATR flights and the ones at which
# VSPaero computes the stability matrices
Velocities=(56.1,71,71,90,130)
rho_vec=(1.225,1.225,1.225,1.225*0.7812,1.225*0.4549)
#a_sound=(340.3,340.3,340.3,330.6,310.2)
Mach=[0.0,0.2,0.3,0.4,0.5]

""" Algorithm start """

#--- List all .stab file from vsp aero and read the coeff ---- 

if g.hangar['aircraft']=='ATR72':
    if g.hangar['version']=='original' or g.hangar['version']=='DEPoriginal':
        """
#        if g.hangar['version']=='original' and g.N_eng!=2:
#            sys.exit("Not possible to use "+str(g.N_eng)+" engines with version "+g.hangar['version'])
        if g.hangar['version']=='original' and g.nofin==True:
            print("WARNING : Using "+g.hangar['version']+" without rudder. Not recommended...")
        if g.hangar['version']=='DEPoriginal' and g.inop!=0 and g.nofin==True:
            print("WARNING : Using "+g.hangar['version']+" with inoperative engine and without rudder. Not recommended...")
        filename=['ATR72_SI_MTOW_Control_flap.stab', 'ATR72_SI_MTOW_Control_flapVmax.stab',
                  'ATR72_SI_MTOW_Control_Vmin.stab','ATR72_SI_MTOW_Control_Manoeuver.stab',
                  'ATR72_SI_MTOW_Control_Cruise.stab']
        Matrix=ReadCoefAero.ReadCoef(filename)
        # --- Seperate the derivatives and the base aero coefficients
        CoefMatrix=Matrix[:,1:]
        BaseAero=Matrix[:,1]
        
#        if g.hangar['version']=='DEPoriginal':
#            g.set_nofin(True)
        """
        if g.hangar['version']=='original' and g.nofin==True:
            print("WARNING : Using "+g.hangar['version']+" without rudder. Not recommended...")
        if g.hangar['version']=='DEPoriginal' and g.inop!=0 and g.nofin==True:
            print("WARNING : Using "+g.hangar['version']+" with inoperative engine and without rudder. Not recommended...")
        filename=['ATR72_SI_MTOW_Control_flap.stab', 'ATR72_SI_MTOW_Control_flapVmax.stab',
                  'ATR72_SI_MTOW_Control_Vmin.stab','ATR72_SI_MTOW_Control_Manoeuver.stab',
                  'ATR72_SI_MTOW_Control_Cruise.stab']
        path='ATR72_SI_MTOW_Control_FinLess_DegenGeom_19_2_15h32/'#'D:/e.nguyen-van/Documents/DEP_control/StabilityMapATR72/ATR72_SI_MTOW_Control_FinLess_DegenGeom_19_2_15h32/'
        filenameNoFin=[path+'_FinLess_Mach0.stab',path+'_FinLess_Mach02.stab',path+'_FinLess_Mach03.stab',
                   path+'_FinLess_Mach04.stab',path+'_FinLess_Mach05.stab']
        MatrixNoFin=ReadCoefAero.ReadCoef(filenameNoFin,col=(11,24,37,50,63,76,128),lign=(47,48,49,50,51,52))
        MatrixFin=ReadCoefAero.ReadCoef(filename)
        # Use coefficient of no fin ATR and add the coefficient for the rudder.
        Matrix=np.copy(MatrixNoFin)
        Matrix[:,-1]=MatrixFin[:,-1]
        #
        print("Adjusting Kf, Kh and VT size")
        print("New Kf and Kh")
        print(g.AdjustVT())
        CoefMatrix=g.NicolosiCoef(Matrix[:,1:], Mach)
        BaseAero=Matrix[:,1]
    
    elif g.hangar['version']=='DEPnofin':
        filename=['ATRfinless_flap.stab', 'ATRfinless_flapVmax.stab',
                  'ATRfinless_Vmin.stab','ATRfinless_Manoeuver.stab','ATRfinless_Cruise.stab']
        Matrix=ReadCoefAero.ReadCoef(filename,col=(11,24,37,50,63,76,128))
        # --- Seperate the derivatives and the base aero coefficients
        CoefMatrix=Matrix[:,1:]
        BaseAero=Matrix[:,1]
        
    
    else:
        sys.exit("Version '"+g.hangar['version']+"' is not a valid case. Stopping")
        

# --- Interpol coefficients for test case ---
# Find sound velocity and air density
atmospher=g.GetAtmo(H_base)
a_sound=atmospher[0]
rho_base=atmospher[1]
M_base=V_base/a_sound
Coef_base=AeroForces.CoefInterpol(M_base, CoefMatrix, Mach)

#--- plot coefficients evolution ---
if plotcoef==True:
    vel=np.linspace(0.0,0.5,30) #in Mach
    Cy_beta=np.array([])
    CL_alpha=np.array([])
    Cy_r=np.array([])
    Cy_p=np.array([])
    Cn_r=np.array([])
    Cl_r=np.array([])
    Cn_beta=np.array([])
    
    for i in range(len(vel)):
        LocalMatrix=AeroForces.CoefInterpol(vel[i], CoefMatrix, Mach)
        Cy_beta=np.append(Cy_beta,LocalMatrix[1,1])
        CL_alpha=np.append(CL_alpha,LocalMatrix[2,0])
        Cy_r=np.append(Cy_r,LocalMatrix[1,4])
        Cy_p=np.append(Cy_p,LocalMatrix[1,2])
        Cn_r=np.append(Cn_r,LocalMatrix[-1,4])
        Cl_r=np.append(Cl_r,LocalMatrix[3,4])
        Cn_beta=np.append(Cn_beta,LocalMatrix[-1,1])
        
    
    plt.figure(1)
    plt.scatter(vel,Cy_beta,s=50,color="0")
    axe1=plt.gca()
    axe1.set_ylim(-0.6,0)
    plt.xlabel('Velocity (Mach)')
    plt.ylabel(r"$Cy_\beta$(/rad)")
    plt.grid()
    plt.tight_layout()
    plt.savefig("CybetaMachChange.pdf")
    
    plt.figure(2)
    plt.scatter(vel,Cy_r,s=50,color="0")
    axe2=plt.gca()
    axe2.set_ylim(0.0,0.65)
    plt.ylabel(r"$Cy_r$ ( /rad)")
    plt.xlabel('Velocity (Mach)')
    plt.grid()
    plt.tight_layout()
    plt.savefig("CyrMachChange.pdf")
    
    plt.figure(3)
    plt.scatter(vel,Cn_r,s=50,color="0")
    axe3=plt.gca()
    axe3.set_ylim(-0.3,0)
    plt.ylabel(r"$Cn_r$ (/rad)")
    plt.xlabel('Velocity (Mach)')
    plt.grid()
    plt.tight_layout()
    plt.savefig("CnrMachChange.pdf")
    
    plt.figure(4)
    plt.scatter(vel,Cn_beta,s=50,color="0")
    axe4=plt.gca()
    axe4.set_ylim(0,0.12)
    plt.ylabel(r"$Cn_\beta$ (/rad)")
    plt.xlabel('Velocity (Mach)')
    plt.grid()
    plt.tight_layout()
    plt.savefig("CnbetaMachChange.pdf")
    
#    plt.figure(4)
#    plt.scatter(vel,CL_alpha,s=100)
#    plt.ylabel(r"$CL_\alpha$ (/°)")
#    plt.xlabel('Velocity (Mach)')
#    plt.grid()

# Initialise test and guess vectors:
if g.nofin==False:
    # x =[alpha, p, q, r, phi, theta, delta_a, delta_e, delta_r, delta_i]
    xtest=np.array([5*math.pi/180, 0,0,0.01, 0.01, 0.01, 0, 0,  0]) # no engines yet
    x0=np.array([5*math.pi/180, 0,0,0, 0.00, 0.00, 0, 0.05, 0])
    bnds=( (-5*math.pi/180,alphamax*math.pi/180), (-0.2,0.2), (-0.2,0.2), (-0.2,0.2), (-phimax/180*math.pi,phimax/180*math.pi), (-30/180*math.pi,30/180*math.pi), (-30/180*math.pi,30/180*math.pi), (-20/180*math.pi,20/180*math.pi), (-deltaRmax/180*math.pi,deltaRmax/180*math.pi))
else:
    # x =[alpha, p, q, r, phi, theta, delta_a, delta_e, delta_i]
    xtest=np.array([5*math.pi/180, 0,0,0.01, 0.01, 0.01, 0, 0]) # no engines yet
    x0=np.array([5*math.pi/180, 0,0,0, 0.00, 0.00, 0, 0.05])
    bnds=( (-5*math.pi/180,alphamax*math.pi/180), (-0.2,0.2), (-0.2,0.2), (-0.2,0.2), (-phimax/180*math.pi,phimax/180*math.pi), (-30/180*math.pi,30/180*math.pi), (-30/180*math.pi,30/180*math.pi), (-20/180*math.pi,20/180*math.pi))

# complete the vectors with engines:
eng_vec=np.array([0.1]*g.N_eng)
xtest=np.append(xtest,eng_vec)
x0=np.append(x0,eng_vec)
bnds_eng=((ThrottleMin,ThrottleMax), (ThrottleMin,ThrottleMax))
for i in range(int(g.N_eng/2)):
    bnds=bnds+bnds_eng

# --- imposed conditions --- 
# fix = [V, beta, gamma, omega, H]
if R==0:
    omega=0
else:
    omega = V_base*math.cos(gamma)/R #omega=Vcos(gamma)/R, 0 if R=0

fixtest=np.array([V_base, beta_base, gamma, omega]) 

# put everything in tuples for passing to functions
diccons=(np.copy(fixtest), np.copy(Coef_base), rho_base, g) #fix, CoefMatrix,Velocities, rho, g
dicfobj=(np.copy(fixtest),rho_base,g)

# --- minimization algorithm ---
k=minimize(e.fobjective, np.copy(x0), args=dicfobj, bounds=bnds, constraints={'type' : 'eq', 'fun' : e.Constraints_DEP, 'args' : diccons}, options={'maxiter':MaxIter}, tol= tolerance)

# print results
print(k)
def printx(x, fix):
    V=fix[0]
    alpha=x[0]/math.pi*180
    beta=fix[1]/math.pi*180
    pqr=x[1:4]/math.pi*180
    phi=x[4]/math.pi*180
    theta=x[5]/math.pi*180
    da=x[6]/math.pi*180
    de=x[7]/math.pi*180
        
    print("\nState vector value:")
    print("V= {0:0.2f}m/s, alpha = {1:0.2f}\xb0, beta={2:0.2f}\xb0, phi={3:0.2f}\xb0, theta={4:0.2f}\xb0".format(V,alpha,beta,phi,theta))
    print("p={0:0.4f}\xb0/s q={1:0.4f}\xb0/s r={2:0.4f}\xb0/s".format(*pqr))
    print("da={0:0.2f}\xb0, de= {1:0.2f}\xb0".format(da,de))
    if g.nofin==False:
        print("dr = {0:0.2f}\xb0".format(x[8]/math.pi*180))

printx(k.x, fixtest)

# check if constraints are validated
constraints_calc=e.Constraints_DEP(k.x,*diccons)
print("\nConstraints")
print(constraints_calc)
# Option plot the thrust distribution and rudder
if HistoGouv==True:
    plt.figure(7)
    plt.subplot(2,1,1)
    plt.bar(g.PosiEng,k.x[-g.N_eng:], color="0.45")
    plt.xlabel("Engine location")
    plt.ylabel("Throttle (-)")
    plt.ylim(0,1)
    plt.xlim(-g.b/2-1,g.b/2+1)
    plt.xticks(np.linspace(-g.b/2,g.b/2,5),('-b/2','-b/4','0','b/4','b/2'))
    plt.grid()
    if g.nofin==False:
        axrud=plt.subplot(2,1,2)
        width=0.5
        axrud.barh(0,k.x[8]/math.pi*180,width,color="0.45")
        axrud.set_yticks([0])
        axrud.set_xlim(-deltaRmax,deltaRmax)
        axrud.set_ylim(-1,1)
        axrud.set_xlabel("Rudder Deflection (°)")
        plt.sca(axrud)
        plt.yticks([0],[r'$\delta_R$'])
        plt.grid()
        plt.tight_layout()
        plt.show()
    else:
        #still show zero deflection
        axrud=plt.subplot(2,1,2)
        width=0.5
        axrud.barh(0,0,width,color="0.45")
        axrud.set_yticks([0])
        axrud.set_xlim(-deltaRmax,deltaRmax)
        axrud.set_ylim(-1,1)
        axrud.set_xlabel("Rudder Deflection (°)")
        plt.sca(axrud)
        plt.yticks([0],[r'$\delta_R$'])
        plt.grid()
        plt.tight_layout()
        plt.show()

# ---- Compute jacobian -----
if gojac==True:
    jac=e.Jac_DEP(k.x,*diccons, 0.01)
    
    LignLat=(1,3,5,6)
    ColLat=(1,3,5,6)
    LignLongi=(0,2,4)
    ColLongi=(0,2,4)
    TransiLat=jac[LignLat,:]
    LatJac=TransiLat[:,ColLat]
    TransiLongi=jac[LignLongi,:]
    LongiJac=TransiLongi[:,ColLongi]
    Lateigvals=scipy.linalg.eigvals(LatJac)
    print("Eigen value :")
    print(Lateigvals)
    # plot eigen values
    plt.figure(5)
    plt.scatter(Lateigvals.real, Lateigvals.imag)
    plt.grid()
    
    #Now check if fin variation is activated
    if goFinVariation == True:
        #set up plot
        eigplot=plt.figure()
        ax = eigplot.add_subplot(111)
        ax.grid()
        MarkerColor=abs(np.linspace(-1,0.1,len(FinRatioVec)))
        # Compute equilibrium and jac at base flight point for different fin size
        for i in range(len(FinRatioVec)):
            #Adjust VT size
            if CstSpan==True and CstA==False:
                #cst span variation
                g.AdjustVTcstSpan(FinRatioVec[i])
            elif CstA==True and CstSpan==False:
                #cst Aspect ratio variation
                g.AdjustVT(FinRatioVec[i])
            else:
                sys.exit("Wrong indication, 'CstA' and 'CstSpan' are either both True or False.")
            
            # Now adjust Matrix
            CoefMatrix=g.NicolosiCoef(Matrix[:,1:], Mach)
            # Interpol coefficients
            Coef_base=AeroForces.CoefInterpol(M_base, CoefMatrix, Mach)
            # re-set optim parameter
            diccons=(np.copy(fixtest), np.copy(Coef_base), rho_base, g) #fix, CoefMatrix,Velocities, rho, g
            # run optim
            k=minimize(e.fobjective, np.copy(x0), args=dicfobj, bounds=bnds, constraints={'type' : 'eq', 'fun' : e.Constraints_DEP, 'args' : diccons}, options={'maxiter':MaxIter}, tol= tolerance)
            if k.success == True:
#                print('Optim succeful')
                #Check constraints again adjust threshold
                check=e.Constraints_DEP(k.x,*diccons)
                err=math.sqrt(np.dot(check,check))
                if err<0.01:
                    #good point compute eigen values
                    jac=e.Jac_DEP(k.x,*diccons, 0.01)
                    TransiLat=jac[LignLat,:]
                    LatJac=TransiLat[:,ColLat]
                    Lateigvals=scipy.linalg.eigvals(LatJac)
                    print(Lateigvals)
                    #plot the thing
                    ax.scatter(Lateigvals.real, Lateigvals.imag, color=str(MarkerColor[i]))
        
#        ax.title('Lateral Eigen Values')

    
# --- display informations --- 
print("\nConfiguration : "+g.hangar['version'])
print("Vertical Tail size ratio : {0:0.2f}".format(g.VTsize))
print("Lateral stability, Cy_beta = {0:0.2f}, Cn_beta = {1:0.3f}, Cn_r = {2:0.2f}".format(CoefMatrix[1,1],CoefMatrix[5,1],CoefMatrix[5,4]))
print("Default conditions : Vel_base = {0:0.1f}m/s, Beta = {1:0.1f}°, gamma={2:0.1f}0, omega={3:0.2f}rad/s, Altitude={4:0.0f}m".format(V_base, beta_base/math.pi*180, gamma/math.pi*180, omega,H_base) )
print("Number of engine : {0} \nNumber of inoperative engines : {1}".format(g.N_eng,g.inop))


""" --- Mapping ---- """
# --- Initialisation by interpreting user commands
underscore=MapName.find("_")    # find the separator "_"
FirstVar=MapName[:underscore]   # get first variable = map ligns
SecondVar=MapName[underscore+1:] # get the second variable = map column


if FirstVar=="Vel" and domap==True:
    MapVariablesData=(MapVmin,MapVmax,MapVstep)
elif FirstVar=="Beta" and domap==True:
    MapVariablesData=(MapBetaMin,MapBetaMax,MapBetaStep)
elif FirstVar=="Gamma" and domap==True:
    MapVariablesData=(MapGammaMin,MapGammaMax,MapGammaStep)
elif FirstVar=="Omega" and domap==True:
    MapVariablesData=(MapOmegaMin,MapOmegaMax,MapOmegaStep)
else:
    print("Invalid lign name, map isn't executed")
    domap=False

if domap==True and SecondVar=="Vel":
    MapVariablesData=MapVariablesData+(MapVmin,MapVmax,MapVstep)
elif domap==True and SecondVar=="Beta":
    MapVariablesData=MapVariablesData+(MapBetaMin,MapBetaMax,MapBetaStep)
elif domap==True and SecondVar=="Gamma":
    MapVariablesData=MapVariablesData+(MapGammaMin,MapGammaMax,MapGammaStep)
elif domap==True and SecondVar=="Omega":
    MapVariablesData=MapVariablesData+(MapOmegaMin,MapOmegaMax,MapOmegaStep)
else:
    print("Invalid column name, map isn't executed")
    domap=False

    
# --- Initialisation completed if successfull start the map ---
if domap==True:
    # Create map variable object
    Var=MapVarClass.MapVar(FirstVar, SecondVar, *MapVariablesData)
    ligns=len(Var.FirstDim)
    col=len(Var.SecondDim)

    #Set up matrix to store flight enveloppe points
    Map=np.zeros( (ligns,col) )
    BparamMatrix=np.zeros( (ligns,col,LenBparam) )
    
    #display more informations
    print("\nRun conditions:")
    print("Sweep : " + FirstVar + " and " + SecondVar)
    
    if gojac==True and storeJac==True:
        dim_eigenvals=7#len(eigvals.real)
        MapJac=np.zeros((ligns,col,dim_eigenvals),dtype=complex)
        print("Go Jacobian")
        
    start_time=time.time()
    for i in range(ligns):
        print("Lign {0} out of {1}".format(i,ligns)) #more info 
        if FirstVar=="Vel":
            # coefficients interpolation is needed if Velocity changes
#            rho_base=AeroForces.InterpolRho(Var.FirstDim[i],rho_vec,Velocities)
            Coef_base=AeroForces.CoefInterpol(Var.FirstDim[i]/a_sound, CoefMatrix, Mach)
        elif FirstVar=="H":
            #update atmosphere settings
            Atmospher=g.GetAtmo(Var.FirstDim[i])
            a_sound=Atmospher[0]
            rho_base=atmospher[1]
            Coef_base=AeroForces.CoefInterpol(Var.FirstDim[i]/a_sound, CoefMatrix, Mach)
        
        for j in range(col):
            if SecondVar=="Vel":
                # coefficients interpolation is needed if Velocity changes
#                rho_base=AeroForces.InterpolRho(Var.SecondDim[j],rho_vec,Velocities)
                Coef_base=AeroForces.CoefInterpol(Var.SecondDim[j]/a_sound, CoefMatrix, Mach)
            elif SecondVar=="H":
                #update atmosphere settings
                Atmospher=g.GetAtmo(Var.SecondDim[i])
                a_sound=Atmospher[0]
                rho_base=atmospher[1]
                Coef_base=AeroForces.CoefInterpol(Var.FirstDim[i]/a_sound, CoefMatrix, Mach)
            
            #Refresh imposed conditions vector
            fix=Var.ReArrangefix(i,j, V_base, beta_base, gamma, omega)
            
            diccons=(np.copy(fix), np.copy(Coef_base), rho_base, g) #fix, CoefMatrix,Mach, rho, g
            dicfobj=(np.copy(fix),rho_base,g)
            
            # minimization algorithm
            k=minimize(e.fobjective, x0, args=dicfobj, bounds=bnds, constraints={'type' : 'eq', 'fun' : e.Constraints_DEP, 'args' : diccons},options={'maxiter':30}, tol= 1e-3 )
            
            if k.success == True:
                #Check constraints again adjust threshold
                check=e.Constraints_DEP(k.x,*diccons)
                err=math.sqrt(np.dot(check,check))
                if err<0.01:
                    Map[i,j]=1
                    #it is a success, for each limiting parameter check whether one is close to limit
                    for l in range(LenSmallBparam):
                        #this is for alpha, phi, da and dR limits
                        iName=BparamName[l]
                        if abs(k.x[BparamXvecPosi[iName]])>BparamLim[iName]*BparamThres[iName]:
                            BparamMatrix[i,j,l]=1
                    for l in range(LenSmallBparam,LenBparam):
                        #this is only for engine limits
                        iName=BparamName[l]
                        if k.x[BparamXvecPosi[iName]]>BparamLim['deltaxmax']*BparamThres['deltaxmax'] :#or k.x[BparamXvecPosi[iName]]<BparamLim['delatxmin']*BparamThres['delatxmin']:
                            BparamMatrix[i,j,l]=1
                    
                    #Compute jac :
                    if gojac==True :
                        jac=e.Jac_DEP(k.x,*diccons,0.01)
                        EigVal=scipy.linalg.eigvals(jac[:,0:7])
                        for l in range(len(EigVal)):
                            if EigVal[l].real>0:
                                Map[i,j]=2
                        if storeJac==True:
                            MapJac[i,j,:]=EigVal

            else:
                Map[i,j]=0
                """#it is a failure, for each limiting parameter check which one is beyond limit
                for l in range(LenSmallBparam):
                    #this is for alpha, phi, da and dR limits
                    iName=BparamName[l]
                    if abs(k.x[BparamXvecPosi[iName]])>BparamLim[iName]*BparamThres[iName]:
                        BparamMatrix[i,j,l]=1
                for l in range(LenSmallBparam,LenBparam):
                    #this is only for engine limits
                    iName=BparamName[l]
                    if k.x[BparamXvecPosi[iName]]>BparamLim['deltaxmax']*BparamThres['deltaxmax'] :#or k.x[BparamXvecPosi[iName]]<BparamLim['delatxmin']*BparamThres['delatxmin']:
                        BparamMatrix[i,j,l]=1
                """
        
    print("--- %s seconds ---" % (time.time()-start_time))
    
    # --- post process the matrix ---
    MarkerStyle=["s","^","v","<",">","+","x","*"]
    lenMarker=len(MarkerStyle)
    LineStyle={'alpha':":",'phi':"-.",'deltaa':"--",'deltaR':"-",'deltax':"--"}
    Labellist={'alpha':"Stall",'phi':"$\phi$>5°",'deltaa':"Aileron saturation",'deltaR':"Rudder saturation",'deltax':"Engine Saturation"}

    Var.Dim2deg()   # transform the variables in degree for plot
    scatter=np.array([]).reshape(0,2)

    # Treat unstable points if necessary
    scatterIm=np.array([0,0])
    for i in range(ligns):
        xscatter=Var.FirstDim[i]
        for j in range(col):
            if Map[i,j]==2:
                scatterIm=np.vstack((scatterIm, np.array([xscatter,Var.SecondDim[j]])))
   
    # Treat bounds
    EngSat={}
    EngSatName=[]
    for i in range(lenMarker):
        EngSatName.append(str(i+1)+"Sat")
        EngSat[EngSatName[i]]=np.array([]).reshape(0,2)

    for i in range(ligns):
        for j in range(col):
            if Map[i,j]==1:
                #Then pull limiting parameters
                for k in range(LenSmallBparam):
                    # This is for alpha, phi or stall
                    if BparamMatrix[i,j,k]==1:
                        kName=BparamName[k]
                        #save the variables values
                        DispBparam[kName]=np.vstack([DispBparam[kName],[Var.FirstDim[i],Var.SecondDim[j]]])

                #Now check Engine saturation
                count=0 #initialize number of engine saturated
                for k in range(LenSmallBparam, LenBparam):
                    #Now check which engine are saturated and count them
                    if BparamMatrix[i,j,k]==1:
                        count+=1
                
                # count is the total number of saturated engine. Save points in EngSat dictionary
                if count==0:
                    scatter=np.vstack([scatter,[Var.FirstDim[i],Var.SecondDim[j]]])
                for k in range(len(EngSatName)):
                    if (k+1)==count:
                        EngSat[EngSatName[k]]=np.vstack([EngSat[EngSatName[k]],[Var.FirstDim[i],Var.SecondDim[j]]])
                    

    #plot the figure
    plt.figure(6)
    ax=plt.subplot()
    xaxes=Var.getXlocator()
    yaxes=Var.getYlocator()
    ax.xaxis.set_minor_locator(ticker.FixedLocator(xaxes))
    ax.yaxis.set_minor_locator(ticker.FixedLocator(yaxes))

    #if beta is one variable, plot fading marker
    if Var.FirstVar=="Beta" or Var.SecondVar=="Beta":
        MarkerColor=["0.45","0.75"]
        if Var.FirstVar=="Beta":
            for l in range(len(scatter[:,0])):
                if scatter[l,0]>-BetaStall and scatter[l,0]<BetaStall:
                    ax.scatter(scatter[l,0],scatter[l,1], color="0.45",s=20)
                else :
                    ax.scatter(scatter[l,0],scatter[l,1], color="0.75",s=20)
            for i in range(len(EngSatName)):
                #Take care of the engine saturation. For each success with at least one engine saturated
                iName=EngSatName[i]
                for l in range(len(EngSat[iName][:,0])):
                    if EngSat[iName][l,0]>-BetaStall and EngSat[iName][l,0]<BetaStall:
                        ax.scatter(EngSat[iName][l,0],EngSat[iName][l,1],marker=MarkerStyle[i], color="0.45")
                    else:
                        ax.scatter(EngSat[iName][l,0],EngSat[iName][l,1],marker=MarkerStyle[i], color="0.75")
        else:
            for l in range(len(scatter[:,1])):
                if scatter[l,1]>-BetaStall and scatter[l,1]<BetaStall:
                    ax.scatter(scatter[l,0],scatter[l,1], color="0.45",s=20)
                else :
                    ax.scatter(scatter[l,0],scatter[l,1], color="0.75",s=20)
            for i in range(len(EngSatName)):
                #Take care of the engine saturation. For each success with at least one engine saturated
                iName=EngSatName[i]
                for l in range(len(EngSat[iName][:,1])):
                    if EngSat[iName][l,1]>-BetaStall and EngSat[iName][l,1]<BetaStall:
                        ax.scatter(EngSat[iName][l,0],EngSat[iName][l,1],marker=MarkerStyle[i], color="0.45")
                    else:
                        ax.scatter(EngSat[iName][l,0],EngSat[iName][l,1],marker=MarkerStyle[i], color="0.75")
    else:
        ax.scatter(scatter[1:,0],scatter[1:,1],color="0.45",s=20)
        for i in range(len(EngSatName)):
            #Take care of the engine saturation. 
            iName=EngSatName[i]
            ax.scatter(EngSat[iName][:,0],EngSat[iName][:,1],marker=MarkerStyle[i], color="0.45")

    if len(scatterIm)>2:
        ax.scatter(scatterIm[1:,0], scatterIm[1:,1],color="r")
    for i in range(LenSmallBparam):
        iName=BparamName[i]
        if len(DispBparam[iName][:,0])>1:
              ax.plot(DispBparam[iName][:,0],DispBparam[iName][:,1],linestyle=LineStyle[iName],color="0",label=Labellist[iName])
    # for i in range(len(EngSatName)):#LenSmallBparam,LenBparam):
    #     #Take care of the engine saturation. For each success with at least one engine saturated
    #     #iName=BparamName[i]
    #     #ax.scatter(DispBparam[iName][:,0],DispBparam[iName][:,1],marker="s",color="0.45")
    #     iName=EngSatName[i]
    #     ax.scatter(EngSat[iName][:,0],EngSat[iName][:,1],marker=MarkerStyle[i], color="0.45")
#    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='0.8', linestyle='-')
    plt.xlabel(Var.getName('x') + " (" +Var.getUnit('x')+")")
    plt.ylabel(Var.getName('y') + " (" +Var.getUnit('y')+")")
    axe=plt.gca()
    axe.set_xlim(Var.FirstDim[0],Var.FirstDim[-1])
    axe.set_ylim(Var.SecondDim[0],Var.SecondDim[-1])
    ax.legend()
    plt.show()
    #plt.savefig(g.hangar['version']+"Map"+MapName+"fin"+str(g.VTsize)+"Eng"+str(g.N_eng+g.inop)+"Rud"+str(g.nofin)+".pdf")


    #plot Jac
    if gojac==True and storeJac==True:
        plt.figure(8)
        plt.grid()
        
        for i in range(ligns):
            eigs=MapJac[i,5,:]
            plt.scatter(eigs.real,eigs.imag)
