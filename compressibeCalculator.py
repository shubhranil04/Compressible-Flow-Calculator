import numpy as np 
from scipy.optimize import fsolve, minimize  
from tabulate import tabulate

'''
    This script calculates compressible aerodynamic relations for 1-D and quasi 1-D 
    isentropic flow, normal shock waves, and oblique shock waves using equations
    from compressible flow theory.
    
    Solution Strategy - 
        (quasi)1-D Isentropic Flow Relations - Provided any of the flow property ratios 
        or angles, the local Mach number is first calculated through iterative methods
        and the remaining property ratios are then computed directly using formulae.
        
        Normal Shock Relations - Provided and of the flow property ratios, the upstream 
        Mach number is first calculated through iterative methods
        and the remaining property ratios are then computed directly using formulae.
        
        Oblique Shock Relations - Provided the upstream Mach number and turn angle,
        shock angle or incident normal Mach number, the shock angle is first calculated 
        through iterative methods   and the remaining property ratios are then 
        computed directly using formulae.
       
    Features:
    	1) User-friendly command line interface for user input.
    	2) Supports several property ratios as inputs.
    	3) Ensures inputs are withing acceptable ranges.
    	
    Dependencies : Numpy, Scipy, Tabulate
        To install the dependencies, use : pip install numpy scipy tabulate
'''

# Isentropic Flow Relations

def ImachAngle (M):
    if M>=1:
        return np.arcsin(1/M)*180/np.pi
    else:
        return None

def IprandtlMeyerAngle (M,gamma): 
    if M>=1:
        return (np.sqrt((gamma+1)/(gamma-1))*np.arctan(np.sqrt((gamma-1)/(gamma+1)*(M**2-1))) - np.arctan(np.sqrt(M**2-1)))*180/np.pi
    else:
        return None
    
def ITbyT0 (M,gamma):
    return 1/(1 + (gamma-1)/2*M**2)

def IPbyP0 (M,gamma):
    return np.power(1 + (gamma-1)/2*M**2,-gamma/(gamma-1))

def IRhobyRho0 (M,gamma):
    return np.power(1 + (gamma-1)/2*M**2,-gamma/(gamma-1))*(1 + (gamma-1)/2*M**2)

def IAbyAcrit (M,gamma):
    return 1/M*np.power(2/(gamma+1)*(1+(gamma-1)/2*M**2),(gamma+1)/2/(gamma-1))

def IPbyPcrit (M,gamma):
    return IPbyP0(M,gamma)/IPbyP0(1,gamma)

def IRhobyRhocrit (M,gamma):
    return IRhobyRho0(M,gamma)/IRhobyRho0(1,gamma)

def ITbyTcrit (M,gamma):
    return ITbyT0(M,gamma)/ITbyT0(1,gamma)

# Normal Shock Relations

def NSM2 (M1,gamma):
    return np.sqrt((1+(gamma-1)/2*M1**2)/(gamma*M1**2-(gamma-1)/2))

def NSRho2byRho1 (M1,gamma):
    return (gamma+1)*M1**2/(2+(gamma-1)*M1**2)

def NSP2byP1 (M1,gamma):
    return 1+2*gamma/(gamma+1)*(M1**2-1)

def NST2byT1(M1, gamma):
    return NSP2byP1(M1,gamma)/NSRho2byRho1(M1,gamma)

def NSP02byP01 (M1,gamma):
    return np.power(((gamma+1)/2*M1**2/(1+(gamma-1)/2*M1**2)),gamma/(gamma-1))/np.power(2*gamma/(gamma+1)*M1**2-(gamma-1)/(gamma+1),1/(gamma-1))

def NSP1byP02 (M1,gamma):
    return IPbyP0(M1,gamma)/NSP02byP01(M1,gamma)

# Oblique shock relations

def OSdelta(M1,gamma,beta):
    beta = np.pi*beta/180  
    return 180/np.pi*np.arctan((2/np.tan(beta)*(M1**2 * (np.sin(beta))**2 - 1))/(2+M1**2*(gamma+np.cos(2*beta))))

def OSminbeta(M1,gamma):
    return fsolve(lambda beta: OSdelta(M1, gamma, beta), 20)

def OSmaxdelta(M1,gamma):
    result = minimize(lambda beta: -OSdelta(M1,gamma,beta), 40)
    maxdelta_beta = result.x[0]
    maxdelta = -result.fun
    return maxdelta, maxdelta_beta

def OSbetaWeak(M1,gamma,delta):
    return fsolve(lambda beta: delta - OSdelta(M1,gamma,beta),20)

def OSbetaStrong(M1,gamma,delta):
    return fsolve(lambda beta: delta - OSdelta(M1,gamma,beta),90)

def OSM1n(M1,gamma,beta):
    return M1*np.sin(beta*np.pi/180)

def OSM2n(M1,gamma,beta):
    return NSM2(OSM1n(M1,gamma,beta),gamma)

def OSM2(M1,gamma,beta):
    return OSM2n(M1,gamma,beta)/np.sin(np.pi/180*(beta-OSdelta(M1,gamma,beta)))

def OSP2byP1(M1,gamma,beta):
    return NSP2byP1(OSM1n(M1,gamma,beta),gamma)

def OSRho2byRho1(M1,gamma,beta):
    return NSRho2byRho1(OSM1n(M1,gamma,beta),gamma)

def OST2byT1(M1,gamma,beta):
    return NST2byT1(OSM1n(M1,gamma,beta),gamma)

def OSP02byP01(M1,gamma,beta):
    M2 = OSM2(M1,gamma,beta)
    return 1/IPbyP0(M2,gamma) * OSP2byP1(M1,gamma,beta) * IPbyP0(M1,gamma)

# Main

print("Welcome to Shubhranil's Compressible Aerodynamics Calculator.")
UserInput = 0

while (UserInput != 4):
    print("\nEnter '1' for Isentropic Flow Relations, '2' for Normal Shock Relations, '3' for Oblique Shock Relations, '4' to Exit the program. ")
    
    UserInput = int(input())
    
    if UserInput == 1:
        print("\nEnter heat capacity ratio (gamma).")
        
        gamma = float(input())
        
        if gamma <= 1:
            print("\nInvalid input. Gamma must be greater than 1.")
            continue
        
        print("\nSelect input : 1-Mach Number, 2-T/T0, 3-P/P0, 4-Rho/Rho0, 5-A/A*(subsonic), 6-A/A*(supersonic), 7-Mach Angle (deg.), 8-Prandtl-Meyer Angle (deg.). ")
        
        inputType = int(input())
            
        if inputType == 1:
            
            print("\nEnter Mach number.")
            
            M = float(input())
            
            if M<=0:
                print("\nInvalid input. Mach number must be greater than 0.")
                
        elif inputType == 2:
            
            print('\nEnter T/T0.')   
            
            TbyT0 = float(input()) 
            
            if (TbyT0 >= 1) or (TbyT0 <= 0):
                print("\nInvalid input. T/T0 must be between 0 and 1.")
                continue
            
            M = fsolve(lambda M: TbyT0 - ITbyT0(M,gamma), 1)
            
        elif inputType == 3:
            
            print('\nEnter P/P0.')   
            
            PbyP0 = float(input()) 
            
            if (PbyP0 >= 1) or (PbyP0 <= 0):
                print("\nInvalid input. P/P0 must be between 0 and 1.")
                continue
            
            M = fsolve(lambda M: PbyP0 - IPbyP0(M,gamma), 1)

        elif inputType == 4:
            
            print('\nEnter Rho/Rho0.')   
            
            RhobyRho0 = float(input()) 
            
            if (RhobyRho0 >= 1) or (RhobyRho0 <= 0):
                print("\nInvalid input. Rho/Rho0 must be between 0 and 1.")
                continue
            
            M = fsolve(lambda M: RhobyRho0 - IRhobyRho0(M,gamma), 1)
            
        elif inputType == 5:
            
            print('\nEnter A/A* (subsonic).')   
            
            AbyAcrit = float(input()) 
            
            if AbyAcrit < 1:
                print("\nInvalid input. A/A* must be greater than 1.")
                continue
            
            M = fsolve(lambda M: AbyAcrit - IAbyAcrit(M,gamma), 0.01)
            
        elif inputType == 6:
            
            print('\nEnter A/A* (supersonic).')   
            
            AbyAcrit = float(input()) 
            
            if AbyAcrit < 1:
                print("\nInvalid input. A/A* must be greater than 1.")
                continue
            
            M = fsolve(lambda M: AbyAcrit - IAbyAcrit(M,gamma), 2)

        elif inputType == 7:
            
            print('\nEnter Mach Angle (deg.).')   
            
            machAngle = float(input()) 
            
            if (machAngle > 90) or (machAngle < 0):
                print("\nInvalid input. Mach Angle must be between 0 and 90 degrees.")
                continue
            
            M = fsolve(lambda M: machAngle - ImachAngle(M), 1)
            
        elif inputType == 8:
            
            print('\nEnter Prandtl-Meyer Angle (deg.).')   
            
            prandtlMeyerAngle = float(input()) 
            
            maxPMAngle = IprandtlMeyerAngle(1e10,gamma)
            
            if (prandtlMeyerAngle < 0) or (prandtlMeyerAngle > maxPMAngle):
                print("Invalid input. Prandtl-Meyer Angle must be between 0 and {:.2f} degrees.".format(maxPMAngle))
                continue
            
            M = fsolve(lambda M: prandtlMeyerAngle - IprandtlMeyerAngle(M,gamma), 1)
        
        else :
            print('\nInvalid input.')
            continue
            
        machAngle = ImachAngle(M)
            
        prandtlMeyerAngle = IprandtlMeyerAngle(M,gamma)
            
        PbyP0 = IPbyP0(M,gamma)
            
        RhobyRho0 = IRhobyRho0(M,gamma)
            
        TbyT0 = ITbyT0(M,gamma)
            
        AbyAcrit = IAbyAcrit(M,gamma)
            
        PbyPcrit = IPbyPcrit(M,gamma)
            
        RhobyRhocrit = IRhobyRhocrit(M,gamma)
            
        TbyTcrit = ITbyTcrit(M,gamma)
        
        headers_1 = ["Mach Number", "Mach Angle", "Prandtl-Meyer Angle", "P/P0", "Rho/Rho0"]
        
        headers_2 = ["T/T0", "P/P*", "Rho/Rho*", "T/T*", "A/A*"]
        
        result_1 = [[M, machAngle, prandtlMeyerAngle, PbyP0, RhobyRho0]]
            
        result_2 = [[TbyT0, PbyPcrit, RhobyRhocrit, TbyTcrit, AbyAcrit]]
            
        print("\n\t\t\t Isentropic Flow Relations \n")
        print(tabulate(result_1, headers = headers_1))
        print("")
        print(tabulate(result_2, headers = headers_2))
    
    elif UserInput == 2:
        print("\nEnter heat capacity ratio (gamma).")
        
        gamma = float(input())
        
        if gamma <= 1:
            print("\nInvalid input. Gamma must be greater than 1.")
            continue
        
        print("\nSelect input : 1-M1, 2-M2, 3-P2/P1, 4-Rho2/Rho1, 5-T2/T1, 6-P02/P01, 7-P1/P02. ")
        
        inputType = int(input())
        
        if inputType == 1:
            print("\nEnter M1.")
            
            M1 = float(input())
            
            if M1 < 1:
                print("\nInvalid input. M1 must be greater than 1.")
                continue
            
        elif inputType == 2:
            print("\nEnter M2.")
            
            M2 = float(input())
            
            minM2 = NSM2(1e10,gamma)
            
            if (M2 > 1) or (M2 < minM2):
                print('\nInvalid input. M2 must be between {:.2f} and 1.'.format(minM2))
                continue
            
            M1 = fsolve(lambda M1: M2 - NSM2(M1,gamma), 2)
            
        elif inputType == 3:
            print("\nEnter P2/P1.")
            
            P2byP1 = float(input())
            
            if P2byP1 < 1:
                print("\nInvalid input. P2/P1 must be greater than 1.")
                continue
            
            M1 = fsolve(lambda M1: P2byP1 - NSP2byP1(M1,gamma),2)
            
        elif inputType == 4:
            print("\nEnter Rho2/Rho1.")
            
            Rho2byRho1 = float(input())
            
            maxRho2byRho1 = NSRho2byRho1(1e10,gamma)
            
            if (Rho2byRho1 < 1) or (Rho2byRho1 > maxRho2byRho1):
                print("\nInvalid input. Rho2/Rho1 must be between 1 and {:.2f}.".format(maxRho2byRho1))
                continue
            
            M1 = fsolve(lambda M1: Rho2byRho1 - NSRho2byRho1(M1,gamma),2)
        
        elif inputType == 5:
            print("\nEnter T2/T1.")
            
            T2byT1 = float(input())
            
            if T2byT1 < 1:
                print("\nInvalid input. T2/T1 must be greater than 1.")
                continue
            
            M1 = fsolve(lambda M1: T2byT1 - NST2byT1(M1,gamma),2)
            
        elif inputType == 6:
            print("\nEnter P02/P01.")
            
            P02byP01 = float(input())
            
            if (P02byP01 >= 1) or (P02byP01 <= 0):
                print("\nInvalid input. P02/P01 must be between 0 and 1.")
                continue
            
            M1 = fsolve(lambda M1: P02byP01 - NSP02byP01(M1,gamma),2)

        elif inputType == 7:
            print("\nEnter P1/P02.")
            
            P1byP02 = float(input())
            
            maxP1byP02 = NSP1byP02(1,gamma)

            if (P1byP02 >= maxP1byP02) or (P1byP02 <= 0):
                print("\nInvalid input. P1byP02 must be between 0 and {:.3f}.".format(maxP1byP02))
                continue
            
            M1 = fsolve(lambda M1: P1byP02 - NSP1byP02(M1,gamma),2)
        
        else :
            print("\nInvalid input.")
            continue
            
        M2 = NSM2(M1,gamma)
        
        P2byP1 = NSP2byP1(M1,gamma)
        
        Rho2byRho1 = NSRho2byRho1(M1,gamma)
        
        T2byT1 = NST2byT1(M1,gamma)
        
        P02byP01 = NSP02byP01(M1,gamma)
        
        P1byP02 = NSP1byP02(M1,gamma)
        
        headers_1 = ["M1", "M2", "P2/P1", "Rho2/Rho1","T2/T1", "P02/P01", "P1/P02"]
        
        result_1 = [[M1, M2, P2byP1, Rho2byRho1,T2byT1, P02byP01, P1byP02]]
            
        print("\n\t\t\t Normal Shock Relations \n")
        print(tabulate(result_1, headers = headers_1))

    elif UserInput == 3:
        print("\nEnter heat capacity ratio (gamma).")
        
        gamma = float(input())
        
        if gamma <= 1:
            print("\nInvalid input. Gamma must be greater than 1.")
            continue
        
        print("\nEnter M1.")
        
        M1 = float(input())
            
        if M1 <= 1:
            print("\nInvalid input. M1 must be greater than 1.")
            continue
        
        deltaMax = OSmaxdelta(M1,gamma)[0]
        deltaMax_beta = OSmaxdelta(M1,gamma)[1]
        
        betaMin = OSminbeta(M1,gamma)[0]
        
        print("\nSelect input : 1-Turn Angle (weak shock), 2-Turn Angle (strong shock), 3-Shock Angle, 4-M1n.")
        
        inputType = int(input())
        
        if inputType == 1:
            print("\nEnter turn angle (weak shock).")
            
            delta = float(input())
            
            if (delta < 0):
                print("\nInvalid input. Turn angle must be greater than 0.")
                continue
            
            if (delta > deltaMax):
                print("\nTurn angle greater than {:.2f}. Detached bow shock appears.".format(deltaMax))
                continue
            
            beta = OSbetaWeak(M1,gamma,delta)
        
        elif inputType == 2:
            print("\nEnter turn angle (strong shock).")
            
            delta = float(input())
            
            if (delta < 0):
                print("\nInvalid input. Turn angle must be greater than 0.")
                continue
            
            if (delta > deltaMax):
                print("\nTurn angle greater than {:.2f}. Detached bow shock appears.".format(deltaMax))
                continue
            
            beta = OSbetaStrong(M1,gamma,delta)
            
        elif inputType == 3:
            print("\nEnter shock angle.")
            
            beta = float(input())

            if beta < betaMin:
                print("\nInvalid input. Shock angle must be greater than Mach angle ({:.2f}).".format(betaMin))
                continue
            
        elif inputType == 4:
            print("\nEnter M1n.")
            
            M1n = float(input())
            
            beta = fsolve(lambda beta: M1n - OSM1n(M1,gamma,beta),40)
        
        else :
            print("\nInvalid input.")
            continue
        
        M2 = OSM2(M1,gamma,beta)
        
        delta = OSdelta(M1,gamma,beta)
        
        P2byP1 = OSP2byP1(M1,gamma,beta)
        
        Rho2byRho1 = OSRho2byRho1(M1,gamma,beta)
        
        T2byT1 = OST2byT1(M1,gamma,beta)
        
        P02byP01 = OSP02byP01(M1,gamma,beta)
        
        M1n = OSM1n(M1,gamma,beta)
        
        M2n = OSM2n(M1,gamma,beta)
        
        headers_1 = ["M1", "M2", "Turn Angle", "Shock Angle", "P2/P1"]
        
        headers_2 = ["Rho2/Rho1", "T2/T1", "P02/P01", "M1n", "M2n"]
        
        result_1 = [[M1, M2, delta, beta, P2byP1]]
        
        result_2 = [[Rho2byRho1, T2byT1, P02byP01, M1n, M2n]]
            
        print("\n\t\t\t Oblique Shock Relations \n")
        print(tabulate(result_1, headers = headers_1))
        print("")
        print(tabulate(result_2, headers = headers_2))
        
    elif UserInput !=4 : 
        print("\nInvalid input.")
