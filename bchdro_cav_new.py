import numpy as np
import copy

# converted from original matlab code 

def bchdro_cav_new(mag, Rrup, Vs30, ZTOR, mechanism, region):

    # input : M, R,vs30,Ztor,mechanism,region----- all  column vector
    # mechanism must be a column of 0 or 1
    # mechansim = 1 for intraslab

    if mechanism == "interface":
        Fs = 0
    elif mechanism == "intraslab":
        Fs = 1

    Vlin = 400
    n = 1.18
    b = -1.186
    c = 1.88
    C4 = 10
    # Vs30_new = min(Vs30,1000)
    Vs30_new = copy.deepcopy(Vs30)
    vssmall = Vs30 <= Vlin
    vsbig = Vs30 > Vlin
    Zsmall = ZTOR <= 100
    Zbig = 1 - Zsmall

    if region == "1_Alaska":
        MB_slab = 7.5
        MB_face = 7.8
        reg = 1
    elif region == "2_Cascadia":
        MB_slab = 7.2
        MB_face = 8.2
        reg = 2
    elif region == "3_CentralAmerica&Mexico":
        MB_slab = 7.5
        MB_face = 7.8
        reg = 3
    elif region == "4_Japan":
        MB_slab = 7.5
        MB_face = 7.8
        reg = 4
    elif region == "5_NewZealand":
        MB_slab = 7.5
        MB_face = 7.8
        reg = 5
    elif region == "6_SouthAmerica":
        MB_slab = 7.5
        MB_face = 7.8
        reg = 6
    elif region == "7_Taiwan":
        MB_slab = 7.5
        MB_face = 7.8
        reg = 7
    elif region == "0_global":
        MB_slab = 7.5
        MB_face = 7.8
        reg = 8

    C1 = Fs * MB_slab + (1 - Fs) * MB_face
    C1slab_C1face = MB_slab - MB_face

    ## regionalized coefficients

    a1 = 5.683074991752335
    a2 = -1.166611231113404
    a3 = 0.02654321348219018
    a4 = 0.8693449770167709
    a5 = 0.5886472382336402
    a6 = -0.0016714118261620763
    a9 = 0.08546533256839378
    a10 = 1.1388432991377258
    a11 = 0.005418097804478372
    a12 = 0.766775749186763
    a13 = -0.07090743422359473
    a14 = -0.11661127502264568

    phi1 = 0.5968277386018089
    tau1 = 0.38099561590236064
    cav = []

    for j, rrup in enumerate(Rrup):

        ## prediction

        Msmall = mag <= C1
        Mbig = mag > C1
        intercept = a1 + a4 * C1slab_C1face * Fs

        f_geom = (a2 + a14 * Fs + a3 * (mag - 7.8)) * np.log(
            rrup + C4 * np.exp((mag - 6) * a9)
        )

        f_attn = a6 * rrup + a10 * Fs

        f_mag = Msmall * (a4 * (mag - C1) + a13 * (10 - mag) ** 2) + Mbig * (
            a5 * (mag - C1) + a13 * (10 - mag) ** 2
        )

        f_depth = Zsmall * a11 * (ZTOR - 60) * Fs + Zbig * a11 * (100 - 60) * Fs

        siterock = (a12 + b * n) * np.log(1100 / Vlin)
        IMrock = np.exp(intercept + f_mag + f_geom + f_depth + f_attn + siterock)

        f_site = vssmall * (
            a12 * np.log(Vs30_new / Vlin)
            - b * np.log(IMrock + c)
            + b * np.log(IMrock + c * (Vs30_new / Vlin) ** n)
        ) + vsbig * (a12 + b * n) * np.log(Vs30_new / Vlin)

        lncav = intercept + f_mag + f_geom + f_depth + f_site + f_attn

        cav.append(np.exp(lncav) * 9.81)
        # cav.append(np.exp(lncav))
        # cav[j] = np.exp(lncav) * 9.81
        # tau(j,i) = tau1
        # phi(j,i) = phi1
        # sigs(j,i) = sqrt(tau(j,i)^2 + phi(j,i)^2);
    return cav
