class DirectFeedback(ProportionalIntegralLoop):
    """
    Direct RF FB (DFB) via generator current (ig) using PI controller [1,2].

    The DirectFeedback inherits from ProportionalIntegralLoop so all
    mandatory parameters from ProportionalIntegralLoop should be passed as
    kwargs when initializing a DirectFeedback object.

    To avoid cavity-beam unmatching (large synchrotron oscilation of beam),
    the CavityResonator generator parameters should be set before
    initialization.

    The loop must be added to the CavityResonator object using:
        cav_res.feedback.append(loop)

    Parameters
    ----------
    DFB_gain : float
        Gain of the DFB.
    DFB_phase_shift : float
        Phase shift of the DFB.
    DFB_sample_num : int, optional
        Sample number to monitor Vc for the DFB.
        The averaged Vc value in DFB_sample_num is monitored.
        Units are in bucket numbers.
        If None, value from the ProportionalIntegralLoop is used.
        Default is None.
    DFB_every : int, optional
        Interval to monitor and change Vc for the DFB.
        Units are in bucket numbers.
        If None, value from the ProportionalIntegralLoop is used.
        Default is None.
    DFB_delay : int, optional
        Loop delay of the DFB.
        Units are in bucket numbers.
        If None, value from the ProportionalIntegralLoop is used.
        Default is None.

    Attributes
    ----------
    phase_shift : float
        Phase shift applied to DFB, defined as psi - DFB_phase_shift.
    DFB_psi: float
        Return detuning angle value with Direct RF feedback in [rad].
    DFB_alpha : float
        Return the amplitude ratio alpha of the DFB.
    DFB_gamma : float
        Return the gain gamma of the DFB.
    DFB_Rs : float
        Return the shunt impedance of the DFB in [ohm].

    Methods
    -------
    DFB_parameter_set(DFB_gain, DFB_phase_shift)
        Set DFB gain and phase shift parameters.
    track()
        Tracking method for the Direct RF feedback.
    DFB_Vg()
        Return the generator voltage main and DFB components in [V].
    DFB_fs()
        Return the modified synchrotron frequency in [Hz].

    References
    ----------
    [1] : Akai, K. (2022). Stability analysis of rf accelerating mode with
    feedback loops under heavy beam loading in SuperKEKB. PRAB, 25(10),
    102002.
    [2] : N. Yamamoto et al. (2023). Stability survey of a double RF system
    with RF feedback loops for bunch lengthening in a low-emittance synchrotron
    ring. In Proc. IPAC'23. doi:10.18429/JACoW-IPAC2023-WEPL161

    """

    def __init__(self,
                 DFB_gain: float,
                 DFB_phase_shift: float,
                 DFB_sample_num: int | None = None,
                 DFB_every: int | None = None,
                 DFB_delay: int | None = None,
                 **kwargs):
        super(DirectFeedback, self).__init__(**kwargs)

        if DFB_delay is not None:
            self.DFB_delay = int(DFB_delay)
        else:
            self.DFB_delay = self.delay

        if DFB_sample_num is not None:
            self.DFB_sample_num = int(DFB_sample_num)
        else:
            self.DFB_sample_num = self.sample_num

        if DFB_every is not None:
            self.DFB_every = int(DFB_every)
        else:
            self.DFB_every = self.every

        record_size = int(np.ceil(self.DFB_delay / self.DFB_every))
        if record_size < 1:
            raise ValueError("Bad parameter set : DFB_delay or DFB_every")

        self.DFB_parameter_set(DFB_gain, DFB_phase_shift)
        if np.sum(np.abs(self.cav_res.beam_phasor)) == 0:
            cavity_phasor = self.cav_res.Vc * np.exp(1j * self.cav_res.theta)
        else:
            cavity_phasor = np.mean(self.cav_res.cavity_phasor_record)
        self.DFB_VcRecord = np.ones(record_size, dtype=complex) * cavity_phasor
        self.DFB_vc_previous = np.ones(self.DFB_sample_num,
                                       dtype=complex) * cavity_phasor

        self.DFB_sample_list = range(0, self.ring.h, self.DFB_every)

    @property
    def DFB_phase_shift(self) -> float:
        """Return DFB phase shift."""
        return self._DFB_phase_shift

    @DFB_phase_shift.setter
    def DFB_phase_shift(self, value: float):
        """Set DFB_phase_shift and phase_shift"""
        self._DFB_phase_shift = value
        self._phase_shift = self.cav_res.psi - value

    @property
    def phase_shift(self) -> float:
        """
        Return phase shift applied to DFB.
        Defined as self.cav_res.psi - self.DFB_phase_shift.
        """
        return self._phase_shift

    @property
    def DFB_psi(self) -> float:
        """
        Return detuning angle value with Direct RF feedback in [rad].

        Fig.4 of ref [1].
        """
        return (np.angle(np.mean(self.cav_res.cavity_phasor_record)) -
                np.angle(np.mean(self.ig_phasor_record)))

    @property
    def DFB_alpha(self) -> float:
        """
        Return the amplitude ratio alpha of the DFB.

        Near Eq. (13) of [1].
        """
        ig_forward = self.ig_phasor_record - self.DFB_ig_phasor
        fac = np.abs(
            np.mean(self.DFB_ig_phasor) / np.mean(ig_forward))
        return 20 * np.log10(fac)

    @property
    def DFB_gamma(self) -> float:
        """
        Return the gain gamma of the DFB.

        Near Eq. (13) of [1].
        """
        ig_forward = self.ig_phasor_record - self.DFB_ig_phasor
        fac = np.abs(
            np.mean(self.DFB_ig_phasor) / np.mean(ig_forward))
        return fac / (1 - fac)

    @property
    def DFB_Rs(self) -> float:
        """
        Return the shunt impedance of the DFB in [ohm].

        Eq. (15) of [1].
        """
        return self.cav_res.Rs / (1 + self.DFB_gamma * np.cos(self.DFB_psi))

    def DFB_parameter_set(self, DFB_gain: float, DFB_phase_shift: float):
        """
        Set DFB gain and phase shift parameters.

        Parameters
        ----------
        DFB_gain : float
            Gain of the DFB.
        DFB_phase_shift : float
            Phase shift of the DFB.

        """
        self.DFB_gain = DFB_gain
        self.DFB_phase_shift = DFB_phase_shift

        if np.sum(np.abs(self.cav_res.beam_phasor)) == 0:
            vc = np.ones(self.ring.h) * self.cav_res.Vc * np.exp(
                1j * self.cav_res.theta)
        else:
            vc = self.cav_res.cavity_phasor_record
        vg_drf = -self.DFB_gain * vc * np.exp(1j * self.phase_shift)
        self.DFB_ig_phasor = self.Vg2Ig(vg_drf)

        self.ig_phasor = self.ig_phasor_record - self.DFB_ig_phasor
        self.init_FFconst()

    def track(self):
        """
        Tracking method for the Direct RF feedback.

        Returns
        -------
        None.

        """
        super(DirectFeedback, self).track(False)

        vc_list = np.concatenate(
            [self.DFB_vc_previous, self.cav_res.cavity_phasor_record])
        self.DFB_ig_phasor = np.roll(self.DFB_ig_phasor, 1)
        for index in self.DFB_sample_list:
            # 2) updating Ig using last item of the list
            vg_drf = -self.DFB_gain * self.DFB_VcRecord[-1] * np.exp(
                1j * self.phase_shift)
            self.DFB_ig_phasor[index:] = self.Vg2Ig(vg_drf)
            # Shift the record
            self.DFB_VcRecord = np.roll(self.DFB_VcRecord, 1)
            # 1) recording Vc
            mean_vc = np.mean(vc_list[index:self.DFB_sample_num + index])
            self.DFB_VcRecord[0] = mean_vc
        # update sample_list for next turn
        self.DFB_sample_list = range(index + self.DFB_every - self.ring.h,
                                     self.ring.h, self.DFB_every)
        # update vc_previous for next turn
        self.DFB_vc_previous = self.cav_res.cavity_phasor_record[
            -self.DFB_sample_num:]

        self.ig_phasor_record = self.ig_phasor + self.DFB_ig_phasor

        self.Ig2Vg()

    def DFB_Vg(self, vc: float = -1) -> tuple[float, float]:
        """Return the generator voltage main and DFB components in [V]."""
        if vc == -1:
            vc = np.mean(self.cav_res.cavity_phasor_record)
        vg_drf = -self.DFB_gain * vc * np.exp(1j * self.phase_shift)
        vg_main = np.mean(self.cav_res.generator_phasor_record) - vg_drf
        return vg_main, vg_drf

    def DFB_fs(self, vg_main: float = -1, vg_drf: float = -1) -> float:
        """Return the modified synchrotron frequency in [Hz]."""
        vc = np.mean(self.cav_res.cavity_phasor_record)
        if vg_drf == -1:
            vg_drf = -self.DFB_gain * vc * np.exp(1j * self.phase_shift)
        if vg_main == -1:
            vg_main = np.mean(self.cav_res.generator_phasor_record) - vg_drf
        vc_fixed = vg_main / (1.0 + self.DFB_gain * np.exp(1j * self.phase_shift))
        vg_sum = np.abs(vc_fixed) * np.sin(np.angle(vc_fixed))
        omega_s = 0
        if (vg_sum) > 0.0:
            omega_s = np.sqrt(self.ring.ac * self.ring.omega1 * (vg_sum) /
                              self.ring.E0 / self.ring.T0)
        return omega_s / 2 / np.pi
