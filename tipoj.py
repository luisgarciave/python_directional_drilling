#%%##importar las librerias que utilizaremos
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
#%% ### Datos de entrada
EWo = -1148; # [ft]
NSo = -1638.3; # [ft]
TVDo = 8500; # [ft]
KOP = 1981; # [ft]
BUR = 2.5; # [°/100 ft]
#%%
L = 100; # [ft]
#%%
DHo = ((EWo)**2 + (NSo)**2)*0.5
R = (180/(BUR*np.pi))*L; # [ft]
if R < DHo:
    CD = DHo - R; #[ft]
else:
    CD =R - DHo; #[ft]
    OD = TVDo - KOP; #[ft]
    OC = ((CD)**2 + (OD)**2)*0.5 #[ft]
    BC = ((OC)**2 - (R)**2)*0.5; #[ft]
    gamma = np.arctan( (CD/OD))*(180/np.pi); #deg
    phi = np.arccos(R/OC)*(180/np.pi);
    if R < DHo:
        theta = 90 - (phi - gamma);
    else:
        theta = 90 - (phi + gamma);
#%% Calculamos los puntos importantes de la tabla
EOB_TVD = KOP + R*np.sin(theta * (np.pi/180)); #[ft]
EOB_DH= R*(1 - np.cos(theta * (np.pi/180))); #[ft]
EOB_MD = KOP + (theta/BUR)*L; # [ft]
MDo = EOB_MD + BC; #[ft]
#%% Se realiza el conte del número de tubos en
# cada sección del pozo
### Nota:
    # NT = número de tubos
    # NTT = número total de tubos
#%% Sección vertical
NT1 = KOP/L; NTT1 = int( np.ceil(NT1));
## np.ceil() realiza redondeo hacia arriba
#%% Sección de construccción
NT2 = (EOB_MD - KOP)/L; NTT2 = int( np.ceil(NT2));
#%% Sección tangente
NT3 = (BC/L); NTT3 = int( np.ceil(NT3) );
#%% Número de tubos totales en toda la trayectoria
NTT = NTT1 + NTT2 + NTT3;
#%% Tabla de puntos importantes

Tabla = np.array ([ ['   ', 'TVD', 'DH', 'MD'],
                    ['LOC', 0.0, 0.0, 0.0],
                    ['KOP', KOP, 0.0, KOP],
                    ['EOB', EOB_TVD, EOB_DH, EOB_MD],
                    ['OBJ', TVDo, DHo, MDo]]);
print(Tabla);
#%% Creamos los arreglos donde guardaremos la solución
MD = np.zeros(NTT+1);
DH = np.zeros(NTT+1);
TVD = np.zeros(NTT+1);
NS = np.zeros(NTT+1);
EW = np.zeros(NTT+1);
INC = np.zeros(NTT+1);
AZI = np.zeros(NTT+1);
#%% Selección de cuadrante
alpha = abs( np.arctan(EWo/NSo)); #rad

if EWo > 0 and NSo > 0:
    #### Cuadrante I
    CuadEW = 1.0;
    CuadNS = 1.0;
    beta = alpha; # rad
elif (EWo > 0 and NSo < 0):
    #### Cuadrante II
    CuadEW = 1.0;
    CuadNS = -1.0;
    beta = np.pi - alpha; # rad
elif (EWo < 0 and NSo < 0):
    #### Cuadrante III
    CuadEW = -1.0;
    CuadNS = -1.0;
    beta = np.pi + alpha; # rad
else:
    CuadEW = -1.0;
    CuadNS = 1.0;
    beta = 2*np.pi - alpha; # rad
    
theta_rad = theta*(np.pi/180);
beta_grad = beta*(180/np.pi);
#%% Calculamos la sección vertical desde LOC hasta el KOP
for i in range (1,NTT1):
    TVD[i] = L + TVD[i-1];
    MD[i] = L + MD[i-1];
    TVD[NTT1] = KOP; MD[NTT1]= KOP;
#%% Calculamos la sección de construcción desde el KOP al EOB
for i in range (NTT1+1, NTT1 + NTT2):
    MD[i] = L + MD[i-1];
    INC[i] = INC[i-1] + BUR;
    TVD[i] = KOP + R*np.sin(INC[i] * np.pi/180);
    DH[i] = R*(1 - np.cos(INC[i] * np.pi/180));
    AZI[i] = beta_grad;
    EW[i] = CuadEW*DH[i]*np.sin(alpha);
    NS[i] = CuadNS*DH[i]*np.cos(alpha);
    
MD[NTT1 +NTT2] = EOB_MD;
TVD[NTT1 +NTT2] = EOB_TVD;
DH[NTT1 +NTT2] = EOB_DH;
AZI[NTT1 +NTT2] = beta_grad;
EW[NTT1 +NTT2] = CuadEW*DH[i]*np.sin(alpha);
NS[NTT1 +NTT2] = CuadNS*DH[i]*np.cos(alpha);
INC[NTT1 +NTT2] = theta
#%% Calculamos la sección tangencial desde el EOB hasta el OBJ
for i in range (NTT1 + NTT2+1, NTT):
    MD[i] = L + MD[i-1];
    INC[i] = theta;
    TVD[i] = L*np.cos( theta_rad ) + TVD[i-1];
    DH[i] = L*np.sin( theta_rad ) + DH[i-1];
    AZI[i] = beta_grad;
    EW[i] = CuadEW*DH[i]*np.sin(alpha);
    NS[i] = CuadNS*DH[i]*np.cos(alpha);

MD[NTT] = MDo;
TVD[NTT] = TVDo;
DH[NTT] = DHo;
AZI[NTT] = beta_grad;
EW[NTT] = CuadEW*DH[NTT]*np.sin(alpha);
NS[NTT] = CuadNS*DH[NTT]*np.cos(alpha);
INC[NTT] = theta;    

#%% Realizamos las gráficas correspondientes con la trayectoria del pozo

mpl.rcParams[ 'legend.fontsize' ] = 10
fig = plt.figure(0)
ax = fig.add_subplot(projection='3d')
ax.plot(EW,NS,TVD, 'k', label = 'Pozo J', linewidth = 5)
plt.gca().invert_zaxis()
ax.legend()
plt.show()