import numpy as np
import os
import shutil
from random import *

def thermo(M, CtoO, T, P, name, cea_path):


    #solar values

    H0=9.257E-01
    He0=7.344E-02
    C0=2.272E-04
    O0=4.535E-04  #REDUCE THIS BY Mg and Si
    N0=6.259E-05
    S0=1.428E-05
    P0=2.688E-07
    Na0=1.846E-06
    K0=1.185E-07
    Fe0=2.69E-05
    Ti0=7.77E-08
    V0=9.257E-09
    Met0=C0+O0+N0+S0+P0+Na0+K0+Fe0+Ti0+V0
    HHe0=1.-Met0
    HetoH=0.0793345576
    H0=HHe0/(HetoH+1)
    He0=HetoH*H0
    M0=Met0/H0

    ###modified
    H=1./(M0*M+HetoH+1)
    He=HetoH*H
    Met=1-H-He
    alpha=Met/Met0
    C=alpha*C0
    O=alpha*O0
    N=alpha*N0
    S=alpha*S0
    PP=alpha*P0  #phosphorous
    Na=alpha*Na0
    K=alpha*K0
    Fe=alpha*Fe0
    Ti=alpha*Ti0
    V=alpha*V0
    CplusO=C+O
    O=CplusO/(1.+CtoO)
    C=CplusO-O
    #O=O-3.28*3.274E-05*alpha#-3.210E-05*alpha #subtracting Mg and Si atomic abund

    #read in template cea2.inp file and replace template
    #values with new values

    H2Oarr=np.zeros(P.shape)
    CH4arr=np.zeros(P.shape)
    COarr=np.zeros(P.shape)
    CO2arr=np.zeros(P.shape)
    NH3arr=np.zeros(P.shape)
    N2arr=np.zeros(P.shape)
    HCNarr=np.zeros(P.shape)
    H2Sarr=np.zeros(P.shape)
    TiOarr=np.zeros(P.shape)
    VOarr=np.zeros(P.shape)
    FeHarr=np.zeros(P.shape)
    C2H2arr=np.zeros(P.shape)
    PH3arr=np.zeros(P.shape)
    C2H6arr=np.zeros(P.shape)
    H2arr=np.zeros(P.shape)
    Hearr=np.zeros(P.shape)
    Naarr=np.zeros(P.shape)
    Karr=np.zeros(P.shape)
    Harr=np.zeros(P.shape)
    MMWarr=np.zeros(P.shape)

    #creating run file
    run=open('run_'+name+'.inp','w+')
    run.writelines('cea2_'+name)
    run.close()

    #creating executable
    exe=open('run_'+name+'.com','w+')
    exe.writelines(cea_path+' < '+'run_'+name+'.inp')
    exe.close()

    os.system('chmod +x '+'run_'+name+'.com')

    #Get thermo.lib and trans.lib
    shutil.copyfile(os.path.join(os.path.dirname(__file__),
                                 'data/thermo.lib'), 'thermo.lib')
    shutil.copyfile(os.path.join(os.path.dirname(__file__),
                                 'data/trans.lib'), 'trans.lib')

    for i in range(P.shape[0]):
        infile=open(os.path.join(os.path.dirname(__file__),
                    'data/cea2_COHNSHeNaKFeTiVP_1x_solar_template_CtoO_0.55.inp'))
        text=infile.readlines()  #elements from column 28:37
        infile.close()
        text[0]=text[0].replace(text[0][28:37],	'{:9.3e}'.format(C))		#'{:9.3e}'.format(Cs)
        text[1]=text[1].replace(text[1][28:37],	'{:9.3e}'.format(O))
        text[2]=text[2].replace(text[2][28:37],	'{:9.3e}'.format(H))
        text[3]=text[3].replace(text[3][28:37],	'{:9.3e}'.format(N))
        text[4]=text[4].replace(text[4][28:37],	'{:9.3e}'.format(S))
        text[5]=text[5].replace(text[5][28:37],	'{:9.3e}'.format(He))
        text[6]=text[6].replace(text[6][28:37],	'{:9.3e}'.format(Na))
        text[7]=text[7].replace(text[7][28:37],	'{:9.3e}'.format(K))
        text[8]=text[8].replace(text[8][28:37],'{:9.3e}'.format(Fe))
        text[9]=text[9].replace(text[9][28:37],'{:9.3e}'.format(Ti))
        text[10]=text[10].replace(text[10][28:37],'{:9.3e}'.format(V))
        text[11]=text[11].replace(text[11][28:37],'{:9.3e}'.format(PP))
        text[13]=text[13].replace(text[13][20:25],'{:9.3e}'.format(P[i]))
        text[14]=text[14].replace(text[14][14:19],'{:6.1f}'.format(T[i]))

        outfile=open('cea2_'+name+'.inp','w+')
        outfile.writelines(text)
        outfile.close()

        #execute
        os.system('./run_'+name+'.com')
        output=open('cea2_'+name+'.out')
        out=output.readlines()  #elements from column 28:37

        sub=' H2O            \n'
        gas=1E-14
        if sub in out:
            loc=out.index(sub)
            gas=float(out[loc+1])
        H2Oarr[i]=gas

        sub=' CH4            \n'
        gas=1E-14
        if sub in out:
            loc=out.index(sub)
            gas=float(out[loc+1])
        CH4arr[i]=gas

        sub=' *CO            \n'
        gas=1E-14
        if sub in out:
            loc=out.index(sub)
            gas=float(out[loc+1])
        COarr[i]=gas

        sub=' *CO2           \n'
        gas=1E-14
        if sub in out:
            loc=out.index(sub)
            gas=float(out[loc+1])
        CO2arr[i]=gas

        sub=' NH3            \n'
        gas=1E-14
        if sub in out:
            loc=out.index(sub)
            gas=float(out[loc+1])
        NH3arr[i]=gas

        sub=' *N2            \n'
        gas=1E-14
        if sub in out:
            loc=out.index(sub)
            gas=float(out[loc+1])
        N2arr[i]=gas

        sub=' HCN            \n'
        gas=1E-14
        if sub in out:
            loc=out.index(sub)
            gas=float(out[loc+1])
        HCNarr[i]=gas

        sub=' H2S            \n'
        gas=1E-14
        if sub in out:
            loc=out.index(sub)
            gas=float(out[loc+1])
        H2Sarr[i]=gas

        sub=' *TiO           \n'
        gas=1E-14
        if sub in out:
            loc=out.index(sub)
            gas=float(out[loc+1])
        TiOarr[i]=gas

        sub=' *VO            \n'
        gas=1E-14
        if sub in out:
            loc=out.index(sub)
            gas=float(out[loc+1])
        VOarr[i]=gas

        sub=' PH3            \n'
        gas=1E-14
        if sub in out:
            loc=out.index(sub)
            gas=float(out[loc+1])
        PH3arr[i]=gas

        sub=' C2H6           \n'
        gas=1E-14
        if sub in out:
            loc=out.index(sub)
            gas=float(out[loc+1])
        C2H6arr[i]=gas

        sub=' C2H2,acetylene \n'
        gas=1E-14
        if sub in out:
            loc=out.index(sub)
            gas=float(out[loc+1])
        C2H2arr[i]=gas

        sub=' FeH            \n'
        gas=1E-14
        if sub in out:
            loc=out.index(sub)
            gas=float(out[loc+1])
        FeHarr[i]=gas

        sub=' *H2            \n'
        gas=1E-14
        if sub in out:
            loc=out.index(sub)
            gas=float(out[loc+1])
        H2arr[i]=gas

        sub=' *He            \n'
        gas=1E-14
        if sub in out:
            loc=out.index(sub)
            gas=float(out[loc+1])
        Hearr[i]=gas

        sub=' *Na            \n'
        gas=1E-14
        if sub in out:
            loc=out.index(sub)
            gas=float(out[loc+1])
        Naarr[i]=gas

        sub=' *K             \n'
        gas=1E-14
        if sub in out:
            loc=out.index(sub)
            gas=float(out[loc+1])
        Karr[i]=gas

        sub=' *H             \n'
        gas=1E-14
        if sub in out:
            loc=out.index(sub)
            gas=float(out[loc+1])
        Harr[i]=gas

        sub=' M, (1/n)  '
        mmw=2.35
        for line in out:
            if sub in line:
                mmw=float(line[13:])
        MMWarr[i]=mmw

    os.system('rm *'+name+'*')
    os.remove('thermo.lib')
    os.remove('trans.lib')
    return H2Oarr, CH4arr, COarr, CO2arr, NH3arr, N2arr, HCNarr, H2Sarr,PH3arr, C2H2arr, C2H6arr, Naarr, Karr, TiOarr, VOarr, FeHarr, Harr,H2arr, Hearr,MMWarr




