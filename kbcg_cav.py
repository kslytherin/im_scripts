import numpy as np

# converted from original matlab code 


def kbcg_cav(M, Rrup, vs30, ZTOR, mechanism, region):
    #  CAV models using KBCG functional form

    # input : M, R,vs30,Ztor,mechanism,region

    # switch mechanism
    if mechanism == "interface":
        Fs = 0
    elif mechanism == "intraslab":
        Fs = 1
    else:
        breakpoint()
    # end

    # vs30= min(vs30,1000)
    k1 = 400
    k2 = -1.311
    c = 1.88
    n = 1.18

    theta10 = 0
    ZB_if = 30
    ZB_slab = 80
    vssmall = vs30 <= k1
    vsbig = vs30 > k1
    if region == "1_Alaska":
        MB_slab = 7.2
        MB_face = 8.4
        reg = 1
    elif region == "2_Cascadia":
        MB_slab = 6.4
        MB_face = 7.7
        reg = 2
    elif region == "3_CentralAmerica&Mexico":
        MB_slab = 6.7
        MB_face = 7.3
        reg = 3
    elif region == "4_Japan":
        MB_slab = 7.6
        MB_face = 8.1
        reg = 4
    elif region == "5_NewZealand":
        MB_slab = 7.5
        MB_face = 8.3
        reg = 5
    elif region == "6_SouthAmerica":
        MB_slab = 7
        MB_face = 8.5
        reg = 6
    elif region == "7_Taiwan":
        MB_slab = 7.1
        MB_face = 7.1
        reg = 7
    elif region == "0_global":
        MB_slab = 7.6
        MB_face = 7.9
        reg = 100

    theta1_if = 1.963526293
    theta1_slab = 2.8
    theta2_if = -1.254378973
    theta2_slab = -1.244868758
    theta3 = 0.028024475
    theta4_if = 1.420617599
    theta4_slab = 1.563987271
    theta5 = 0.747139186
    theta_62 = -0.001999881
    theta7 = 0.939084461
    theta9_if = 0.025040951
    theta9_slab = 0.003009148
    theta_nft1 = 0.844321045
    theta_nft2 = 0.196390595
    dZB_if = 0.105201075
    dZB_slab = -0.00771224

    cav = []

    for j, rrup in enumerate(Rrup):
        intercept = (1 - Fs) * theta1_if + Fs * theta1_slab

        f_mag = (1 - Fs) * lh(
            M, MB_face, theta4_if * (MB_face - 6), theta4_if, theta5, 0.1
        ) + Fs * lh(M, MB_slab, theta4_slab * (MB_slab - 6), theta4_slab, theta5, 0.1)

        f_geom = (1 - Fs) * (theta2_if + theta3 * M) * np.log(
            rrup + 10 ** (theta_nft1 + theta_nft2 * (M - 6))
        ) + Fs * (theta2_slab + theta3 * M) * np.log(
            rrup + 10 ** (theta_nft1 + theta_nft2 * (M - 6))
        )

        f_depth = (1 - Fs) * lh(
            ZTOR,
            ZB_if + dZB_if,
            theta9_if * (ZB_if + dZB_if - 15),
            theta9_if,
            theta10,
            1,
        ) + Fs * lh(
            ZTOR,
            ZB_slab + dZB_slab,
            theta9_slab * (ZB_slab + dZB_slab - 50),
            theta9_slab,
            theta10,
            1,
        )

        f_attn = theta_62 * rrup

        siterock = (theta7 + k2 * n) * np.log(1100 / k1)
        PGArock = np.exp(intercept + f_mag + f_geom + f_depth + f_attn + siterock)

        f_site = vssmall * (
            theta7 * np.log(vs30 / k1)
            + k2 * np.log(PGArock + c * (vs30 / k1) ** n)
            - k2 * np.log(PGArock + c)
        ) + vsbig * (theta7 + k2 * n) * np.log(vs30 / k1)

        y = intercept + f_mag + f_geom + f_depth + f_site + f_attn
        # breakpoint()
        # cav[j] = np.exp(y) * 9.81
        cav.append(np.exp(y) * 9.81)

        # tau[j,i] = 0.32
        # phi[j,i] = 0.57
        # sigma[j,i] = sqrt(tau(j,i)^2 + phi(j,i)^2)

    return cav


def lh(x, x0, a, b0, b1, delta):
    res = a + b0 * (x - x0) + (b1 - b0) * delta * np.log(1 + np.exp((x - x0) / delta))

    return res
