from numpy import sin, sinh, sqrt, linspace, zeros, savetxt, array

def distance_integral_element(Input_Param,z):
    dl = sqrt( Input_Param.Omega_M*(1.0+z)**3 +
               Input_Param.Omega_DE*(1.0+z)**(3.0*(1.0+Input_Param.w)) +
               Input_Param.Omega_R*(1.0+z)**4 -
               Input_Param.Omega_k*(1.0+z)**2 )
    return 1.0/dl

def Proper_Distance_Tabulate(Input_Param, z_max):

    # --- convert redshift to proper distance
    integral_n = 10000
    z = linspace(0.0, z_max, integral_n+1)
    dz= z[1]-z[0]
    rz= zeros(len(z))

    integral_dl = 0.0
    if (Input_Param.Omega_k > 0.0):
        for i in range(1,integral_n+1):
            integral_dl += dz*( distance_integral_element(Input_Param,z[i]) +
                                distance_integral_element(Input_Param,z[i-1]) ) / 2.0
            rz[i] = 1.0/(sqrt(abs(Input_Param.Omega_k))) * sin(sqrt(abs(Input_Param.Omega_k))*integral_dl)
    elif (Input_Param.Omega_k == 0.0):
        for i in range(1,integral_n+1):
            integral_dl += dz*( distance_integral_element(Input_Param,z[i]) +
                                distance_integral_element(Input_Param,z[i-1]) ) / 2.0
            rz[i] = integral_dl
    else:
        for i in range(1,integral_n+1):
            integral_dl += dz*( distance_integral_element(Input_Param,z[i]) +
                                distance_integral_element(Input_Param,z[i-1]) ) / 2.0
            rz[i] = 1.0/(sqrt(abs(Input_Param.Omega_k))) * sinh(sqrt(abs(Input_Param.Omega_k))*integral_dl)

    rz = rz * 2997.92458

    # saving:
    f = open("./Output/tabulated_data/Proper_Distance.txt", "w")
    f.write("#    redshift             Proper distance (Mpc/h)\n")        # column names
    savetxt(f, array([z, rz]).T)
