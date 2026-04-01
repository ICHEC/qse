def blockade_radius(rabi_frequency, c6=5420158.53):
    r"""
    Calculate the blockade radius based on the Rabi frequency and C6 coefficient.
    The blockade radius is a measure of the distance at which the Rydberg blockade
    effect occurs, preventing simultaneous excitation of nearby atoms.

    Parameters
    ----------
    rabi_frequency : float
        The Rabi frequency (in rads/μs) associated with the atomic transition.
    c6 : float, optional
        The C6 coefficient (in rads/μs * μm^6) for the van der Waals interaction.
        Default is 5420158.53.

    Returns
    -------
    float
        The blockade radius (in μm).

    Notes
    -----
    The blockade radius is calculated as:

    .. math::
        r_b = \left( \frac{C_6}{\Omega} \right)^{1/6}

    where :math:`\Omega` is the Rabi frequency and :math:`C_6` is the van der Waals
    coefficient.
    """
    return (c6 / rabi_frequency) ** (1 / 6)
