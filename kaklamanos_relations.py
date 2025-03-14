"""
Calculates EQ parameters using relations derived in Kaklamanos (2011)
Started Apr 16, 2024
"""

import numpy as np
import mytools


def calc_width(M, fault_type):
    if isinstance(fault_type, str):
        new_fault_type = np.full(np.shape(M), fault_type, str)
    else:
        new_fault_type = fault_type
    W = 10 ** (-0.76 + 0.27 * M)
    W[new_fault_type == "R"] = 10 ** (-1.61 + 0.41 * M[new_fault_type == "R"])
    W[new_fault_type == "N"] = 10 ** (-1.14 + 0.35 * M[new_fault_type == "N"])
    return W


def calc_zhyp(M, fault_type):
    if isinstance(fault_type, str):
        new_fault_type = np.full(np.shape(M), fault_type, str)
    zhyp = 5.63 + 0.68 * M
    zhyp[new_fault_type == "R" or new_fault_type == "N"] = (
        11.24 - 0.2 * M[new_fault_type == "R" or new_fault_type == "N"]
    )
    zhyp[new_fault_type == "IDK"] = 7.08 + 0.61 * M[new_fault_type == "IDK"]
    return zhyp


def calc_ry(rx, rjb, pre_alpha_deg):
    alpha_deg = np.full_like(rjb, pre_alpha_deg)
    alpha = np.full_like(rjb, np.radians(alpha_deg))
    ry = np.abs(rx / np.tan(np.radians(alpha_deg)))
    ry[np.abs(alpha_deg) == 90] = 0
    acon = np.logical_or(np.abs(alpha_deg) == 180, alpha_deg == 0)
    ry[acon] = rjb[acon]
    return ry


def calc_rrup(rx, ry, rjb, pre_dip_deg, pre_W, pre_ztor):
    W = np.full_like(rx, pre_W)
    ztor = np.full_like(rx, pre_ztor)
    dip_deg = np.full_like(rx, pre_dip_deg)
    delta = np.full_like(rx, np.radians(dip_deg))
    rrup = np.sqrt(rjb**2 + ztor**2)
    con1 = rx < ztor * np.tan(delta)
    con2 = rx > ztor * np.tan(delta) + W / np.cos(delta)
    rrup_p = rx * np.sin(delta) + ztor * np.cos(delta)
    rrup_p[con1] = np.sqrt(rx[con1] ** 2 + ztor[con1] ** 2)
    rrup_p[con2] = np.sqrt(
        (rx[con2] - W[con2] * np.cos(delta[con2])) ** 2
        + (ztor[con2] + W[con2] * np.sin(delta[con2])) ** 2
    )
    rrup[dip_deg != 90] = np.sqrt(rrup_p[dip_deg != 90] ** 2 + ry[dip_deg != 90] ** 2)
    return rrup


def calc_rx(rrup, rjb, pre_dip_deg, W, ztor, pre_alpha_deg):
    dip_deg = np.full_like(W, pre_dip_deg)
    alpha_deg = np.full_like(W, pre_alpha_deg)
    delta = np.full_like(W, np.radians(dip_deg))
    alpha = np.full_like(W, np.radians(alpha_deg))
    rx = (
        rjb
        * np.tan(alpha)
        * np.cos(alpha - np.arcsin(W * np.cos(delta) * np.cos(alpha) / rjb))
    )

    con1 = np.logical_and(
        0 <= alpha_deg,
        np.logical_and(
            alpha_deg < 90, rjb * np.abs(np.tan(alpha)) <= W * np.cos(delta)
        ),
    )

    rx[con1] = rjb[con1] * np.abs(np.tan(alpha[con1]))

    a90con1 = np.logical_and(alpha_deg == 90, rjb > 0)
    rx[a90con1] = rjb[a90con1] + W[a90con1] * np.cos(delta[a90con1])

    a90con2 = np.logical_and(
        alpha_deg == 90, np.logical_and(rjb == 0, rrup < (ztor / np.cos(delta)))
    )
    rx[a90con2] = np.sqrt(rrup[a90con2] ** 2 + ztor[a90con2] ** 2)

    a90con3 = np.logical_and(
        alpha_deg == 90, np.logical_and(rjb == 0, rrup >= (ztor / np.cos(delta)))
    )
    rx[a90con3] = rrup[a90con3] / np.sin(delta[a90con3]) - ztor[a90con3] / np.tan(
        delta[a90con3]
    )

    vertF_con = dip_deg == 90
    rx[vertF_con] = rjb[vertF_con] * np.sin(alpha[vertF_con])

    FW_con = np.logical_and(-180 <= alpha_deg, alpha_deg < 0)
    rx[FW_con] = rjb[FW_con] * np.sin(alpha[FW_con])

    rx[rjb == 0] = 0.5 * W[rjb == 0] * np.cos(delta[rjb == 0])

    return rx


def calc_rx2(rrup, rjb, pre_dip_deg, W, ztor, pre_alpha_deg):
    # Warning!!, this always assumes a positive value for rx
    alpha_deg = np.full_like(W, pre_alpha_deg)
    alpha = np.full_like(W, np.radians(alpha_deg))
    dip_deg = np.full_like(W, pre_dip_deg)
    delta = np.full_like(W, np.radians(dip_deg))
    # fault at an angle and assumming pre_alpha_deg==90?
    rx = W * np.cos(delta) + rjb
    # for vertical fault or on otherside of fault
    rxsqared = rrup[dip_deg == 90] ** 2 - ztor[dip_deg == 90] ** 2
    rx[dip_deg == 90] = np.sqrt(rxsqared)
    # weird conditiion needs to be discussed
    rx[rxsqared < 0] = W[rxsqared < 0] * np.cos(delta[rxsqared < 0]) + rjb[rxsqared < 0]
    rx[dip_deg == 0] = rjb[dip_deg == 0]
    return rx


def calc_alpha(rx, ry):
    alpha = np.full_like(ry, np.degrees(np.arctan(rx / ry)))
    alpha[np.logical_and(rx > 0, ry == 0)] = 90
    alpha[np.logical_and(rx < 0, ry == 0)] = -90
    alpha[np.logical_and(rx == 0, ry > 0)] = 0
    alpha[np.logical_and(rx == 0, ry < 0)] = 180
    alpha[np.logical_and(rx < 0, ry < 0)] = alpha[np.logical_and(rx < 0, ry < 0)] - 180
    alpha[np.logical_and(rx > 0, ry < 0)] = alpha[np.logical_and(rx > 0, ry < 0)] + 180
    return alpha
