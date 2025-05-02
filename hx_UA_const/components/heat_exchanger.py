# components/heat_exchanger.py
import numpy as np
from ..core.sim_cycle import SimCycle

class HeatExchanger:
    """
    General heat exchanger supporting both condenser and evaporator behavior,
    with optional charge mass calculation.
    """
    def __init__(self,
                 sim: SimCycle,
                 N: int,
                 UA_elem: float,
                 T_air: float,
                 mode: str,
                 charge_mode: bool = False,
                 V_elem: float = 0.0):
        assert mode in ('cond', 'eva'), "mode must be 'cond' or 'eva'"
        self.sim = sim
        self.N = N
        self.UA = UA_elem
        self.T_air = T_air
        self.mode = mode
        self.charge_mode = charge_mode
        self.V_elem = V_elem if charge_mode else 0.0

    def exchange(self,
                 mdot: float,
                 P: float,
                 h_in: float):
        h_elem = np.empty(self.N + 1)
        s_elem = np.empty(self.N)
        T_elem = np.empty(self.N)
        q_elem = np.empty(self.N)
        Cp_elem = np.empty(self.N)
        if self.charge_mode:
            rho_elem = np.empty(self.N)
        h_elem[0] = h_in

        for i in range(self.N):
            # get temperature and quality (and density if charge_mode)
            if self.charge_mode:
                T_elem[i], s_elem[i], q_elem[i], rho_elem[i] = self.sim.get_multiple(
                    'HP_inputs', h_elem[i], P, ('T', 'S','Q', 'D'))
            else:
                T_elem[i], s_elem[i], q_elem[i] = self.sim.get_multiple(
                    'HP_inputs', h_elem[i], P, ('T', 'S', 'Q'))

            delta_T = abs(T_elem[i] - self.T_air)
            if 0 < q_elem[i] < 1:
                Q = self.UA * delta_T
            else:
                Cp_elem[i] = self.sim.get_single_no_update('C')
                eps = 1 - np.exp(- self.UA / (mdot * Cp_elem[i]))
                Q = eps * mdot * Cp_elem[i] * delta_T

            # update enthalpy based on mode
            if self.mode == 'cond':
                h_elem[i + 1] = h_elem[i] - Q / mdot
            else:
                h_elem[i + 1] = h_elem[i] + Q / mdot

        h_out = h_elem[-1]
        s_out = s_elem[-1]
        T_out = self.sim.get_single('HP_inputs', h_out, P, ('T'))

        if self.charge_mode:
            m_elem = self.V_elem * rho_elem
            m_tot = np.sum(m_elem)
            return h_out, s_out, T_out, m_tot
        return h_out, s_out, T_out

class Condenser(HeatExchanger):
    def __init__(self,
                 sim: SimCycle,
                 N: int,
                 UA_elem: float,
                 T_air: float,
                 V_elem: float = None):
        charge_mode = V_elem is not None
        super().__init__(sim, N, UA_elem, T_air,
                         mode='cond',
                         charge_mode=charge_mode,
                         V_elem=V_elem or 0.0)

    def exchange(self,
                 mdot: float,
                 P: float,
                 h_in: float):
        return super().exchange(mdot, P, h_in)

class Evaporator(HeatExchanger):
    def __init__(self,
                 sim: SimCycle,
                 N: int,
                 UA_elem: float,
                 T_air: float,
                 V_elem: float = None):
        charge_mode = V_elem is not None
        super().__init__(sim, N, UA_elem, T_air,
                         mode='eva',
                         charge_mode=charge_mode,
                         V_elem=V_elem or 0.0)

    def exchange(self,
                 mdot: float,
                 P: float,
                 h_in: float):
        return super().exchange(mdot, P, h_in)
