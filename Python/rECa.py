#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Import dependencies
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import epsilon_0
from scipy import interpolate
from shapely.geometry import LineString
import sys

# Import FDEM1D (http://github.com/dhanssens)
import FDEM1D


def rECa(sensor, QP_data, IP_data, precision=.001, noise=0, ref_ECa=None, ori_MSa=0, alt_MSa=0, max_ECa=4):
    """
        Calculates the rECa (S/m) of an FDEM (QP and IP) dataset (ppm)

        Parameters
        ----------
        sensor: object
            Sensor object (FDEM1D.Sensor)

        QP_data: np.array
            QP (quadrature-phase or out-of-phase) data (ppm)

        IP_data: np.array
            IP (in-phase) data (ppm)

        precision: float, (optional)
            Approximated required ECa precision (S/m), .001 by default

        noise: float, (optional)
            Instrument noise level (ppm), 0 by default

        ref_ECa: float, (optional)
            Additional reference ECa estimation (S/m), None by default such that EG2015 ECa (Appendix B)
            algorithm is used to estimate the additional ECa value

        ori_MSa: float, (optional)
            Homogeneous half-space MS (-) estimation used to generate the original ECa-QP curve (-), 0 by default

        alt_MSa: float, (optional)
            Altered homogeneous half-space MS estimation used to generate the alternative ECa-QP_a_l_t curve, 0 by default

        max_ECa: float, (optional)
            Maximum sampled homogeneous half-space EC value (S/m)

        Returns
        -------
        rECa: np.array
            Robust apparent electrical conductivity (rECa) (S/m)

        is_ECa_robust: np.array, boolean
            Assessment of the robustness of the rECa values

        Cite
        ----
        Hanssens, D., Delefortrie, S., Bobe, C., De Smedt, P., and M. Van Meirvenne, 2018,
        Robust apparent electrical conductivity (rECa) estimation in ground-based frequency
        domain electromagnetics: Submitted to Geoderma, 2018

        :AUTHOR: Daan Hanssens
        :CONTACT: daan.hanssens@ugent.be
        :REQUIRES: numpy, FDEM (http://github.com/dhanssens), shapely, scipy
    """

    # Get extra information about input structure
    shape = np.shape(QP_data)
    size = np.size(QP_data)

    # Flatten initial matrices/arrays
    rQPdata = np.asarray(QP_data).flatten()
    rIPdata = np.asarray(IP_data).flatten()

    # Generate EC-QP curve
    [EC_range, FWD_QP, FWD_IP, non_robust_ECa] = ECa_QP_curve(sensor, precision=precision, noise=noise, max_ECa=max_ECa, ori_MSa=ori_MSa, alt_MSa=alt_MSa)

    # Check monotony
    if np.all(np.diff(FWD_QP) > 0):

        # Interpolation (pchip)
        iECa = interpolate.pchip_interpolate(FWD_QP, EC_range, rQPdata)

    elif np.all(np.diff(FWD_QP) < 0):

        # Interpolation (pchip) of flipped data, pchip requires sorted data
        iECa = interpolate.pchip_interpolate(np.flip(FWD_QP, axis=0), np.flip(EC_range, axis=0), rQPdata)

    else:

        # Case no initial ECa estimation is included
        if ref_ECa is None:

            # Initialize np.array
            ref_ECa = np.zeros(size)

            # Loop data TODO vectorize
            for ii in range(size):

                if np.isnan(rQPdata[ii]):

                    # Assign NaN value
                    ref_ECa[ii] = np.NaN

                else:

                    # Calculate initial ECa
                    E = 1 / 2 * np.sqrt((FWD_QP - rQPdata[ii]) ** 2 + (FWD_IP - rIPdata[ii]) ** 2)
                    ref_ECa[ii] = EC_range[E.argmin()]

        else:

            # Create init ECa vector
            ref_ECa = np.full(size, ref_ECa)

        # Initialize curve shape
        EC_QP_curve = LineString(list(zip(EC_range.tolist(), FWD_QP.tolist())))
        iECa = np.ones(size) * np.NaN


        # Determine rECa (S/m) TODO vectorize
        for ii in range(size):

            # Check if value exists
            if np.isnan(rQPdata[ii]):

                # Set to NaN
                iECa[ii] = np.NaN

            else:

                # Calculate intersection
                iLine = LineString([(EC_range[0], rQPdata[ii]), (EC_range[-1], rQPdata[ii])])
                mItx = np.asarray(EC_QP_curve.intersection(iLine))

                # Case no intersection
                if not mItx.any():

                    # Set to NaN
                    mItx = np.NaN

                else:

                    # Get relevant data
                    if mItx.ndim > 1:

                        # Grab data
                        mItx = mItx[:, 0]

                        # Calculate nearest
                        near = np.abs(mItx - ref_ECa[ii])

                        # Get iECa (S/m), i.e. get nearest
                        iECa[ii] = mItx[near.argmin()]

                    else:

                        # Get iECa (S/m)
                        iECa[ii] = mItx[0]

                # Clear variable
                del mItx

    # Check robustness
    is_ECa_robust = np.logical_not(np.any(np.abs(np.tile(non_robust_ECa, (size, 1)) - np.tile(iECa.reshape((size, 1)), (1, non_robust_ECa.size))) < precision, axis=1))

    # Reshape to original dimensions
    rECa = np.reshape(iECa, shape)

    # Return output
    return rECa, is_ECa_robust


def ECa_QP_curve(sensor, precision=.001, noise=0, max_ECa=4, min_ECa=0.0001, ori_MSa=.0, alt_MSa=0):
    """
        Calculates the ECa-QP and -IP curve.

        Parameters
        ----------
        sensor: object
            Sensor object (FDEM1D.Sensor)

        precision: float, optional
            Approximated required ECa precision (S/m)

        noise: float, optional
            FDEM instrument noise (ppm)

        max_ECa: float, optional
            Maximum sampled homogeneous half-space EC value (S/m)

        min_ECa: float, optional
            Minimum sampled homogeneous half-space EC value (S/m)

        ori_MSa: float, optional
            Homogeneous half-space MS (-) value used to generate the original ECa-QP curve, 0 by default

        alt_MSa: float, optional
            Altered homogeneous MS value used to generate the alternative ECa-QP_a_l_t curve, 0 by default

        Returns
        -------
        ECa: np.array
            ECa (Apparent electrical conductivity) (S/m)

        QP: np.array
            QP (quadrature-phase or out-of-phase) data (ppm)

        IP: np.array
            IP (in-phase) data (ppm)

        non_robust_ECa: np.array
            Non-robust ECa values (S/m)

        Cite
        ----
        Hanssens, D., Delefortrie, S., Bobe, C., De Smedt, P., and M. Van Meirvenne, 2018,
        Robust apparent electrical conductivity (rECa) estimation in frequency domain
        electromagnetics (including MATLAB, Python and pseudo-code): Submitted to
        Computers & Geosciences, 2018

        :AUTHOR: Daan Hanssens
        :CONTACT: daan.hanssens@ugent.be
        :REQUIRES: numpy, FDEM1D (http://github.com/dhanssens), scipy

        TODO fully paralize FWD model
    """

    # Initialize model characteristics
    sus = np.array([ori_MSa])
    perm = np.array([epsilon_0])
    thick = np.array([0])

    # Initial ECa range (S/m)
    number_of_EC_samples = int(np.round((max_ECa - min_ECa) * 1000))
    if number_of_EC_samples > 20000: raise ValueError('max ECa is too high, make sure it is in S/m')  # Prevent overflow
    ECa = np.linspace(min_ECa, max_ECa, number_of_EC_samples)

    # Initialize QP and IP response arrays
    IP = np.zeros(number_of_EC_samples)
    QP = np.zeros(number_of_EC_samples)

    # Loop samples
    for ii in range(number_of_EC_samples):

        # Update homogeneous half-space EC (S/m)
        con = np.array([ECa[ii]])
        model = FDEM.MODEL.Model(thick, sus, con, perm)

        # Calculate forward response (ppm)
        [IP[ii], QP[ii]] = FDEM1D.Calculate(sensor, model).forward()

    # Assess robustness of curve (if noise is present)
    if noise != 0:

        # Transform precision to indices
        precision_ind = int(precision * 1000)

        # Check if magnetic effects should be accounted for
        if alt_MSa != 0:

            # Initialize response arrays
            IP_alt = np.zeros(number_of_EC_samples)
            QP_alt = np.zeros(number_of_EC_samples)

            # Alter MS (-), calculate magnetic effect
            sus = sus + alt_MSa

            # Loop samples
            for ii in range(number_of_EC_samples):

                # Update homogeneous half-space EC (S/m)
                con = np.array([ECa[ii]])
                model = FDEM.MODEL.Model(thick, sus, con, perm)

                # Calculate forward response (ppm)
                [IP_alt[ii], QP_alt[ii]] = FDEM1D.Calculate(sensor, model).forward()

            # Calculate slope for required precision and magnetic effects
            QP_diff = np.abs(np.diff(QP[::precision_ind])) - (np.abs(QP[::precision_ind] - QP_alt[::precision_ind]))[:-1]

        else:

            # Calculate slope for required precision
            QP_diff = np.abs(np.diff(QP[::precision_ind]))

        # Account for FDEM instrument noise
        mask = QP_diff < noise

        # Grab non-robust ECa and QP values
        non_robust_ECa = ECa[:-precision_ind:precision_ind][mask]

    else:

        # Create vector, all values are de facto robust because 0 ppm noise
        non_robust_ECa = np.ones(number_of_EC_samples)

    # Return ECa, QP, IP, non_robust_ECa
    return ECa, QP, IP, non_robust_ECa
