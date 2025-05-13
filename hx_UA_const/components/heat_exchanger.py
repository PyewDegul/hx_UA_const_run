# components/heat_exchanger.py
import numpy as np
from ..core.sim_cycle import SimCycle
""" 새로운 class 만들기, condenser, evaporator, 단상, 2상에 대한 Nu 관계식 대입 with no pressure drop"""

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
                 V_elem: float):
        assert mode in ('cond', 'eva'), "mode must be 'cond' or 'eva'"
        self.sim = sim
        self.N = N
        self.UA = UA_elem
        self.T_air = T_air
        self.mode = mode
        self.V_elem = V_elem

    def exchange(self,
                 mdot: float,
                 P: float,
                 h_in: float):
        h_elem = np.empty(self.N + 1)
        s_elem = np.empty(self.N)
        T_elem = np.empty(self.N)
        q_elem = np.empty(self.N)
        Cp_elem = np.empty(self.N)
        rho_elem = np.empty(self.N)
        h_elem[0] = h_in

        for i in range(self.N):
            # get temperature and quality (and density if charge_mode)
            T_elem[i], s_elem[i], q_elem[i], rho_elem[i] = self.sim.get_multiple(
                'HP_inputs', h_elem[i], P, ('T', 'S', 'Q', 'D'))

            # 온도차 계산
            delta_T = abs(T_elem[i] - self.T_air)
            if 0 < q_elem[i] < 1:
                # Two phase에서 Const UA 가정
                Q = self.UA * delta_T
            else:
                # C_r == 0 에서 eps-NTU 모델 사용(Al exchanger)
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

        if self.V_elem == 0.0:
            # No charge mode
            return h_out, s_out, T_out
        else:
            m_elem = self.V_elem * rho_elem
            m_tot = np.sum(m_elem)
            return h_out, s_out, T_out, m_tot
        

class Condenser(HeatExchanger):
    def __init__(self,
                 sim: SimCycle,
                 params
                 ):
        super().__init__(sim, params.N_cond, params.UA_cond_elem, params.T_cond_air,
                    mode='cond',
                    V_elem=params.V_elem_cond)
        

    def exchange(self,
                 mdot: float,
                 P: float,
                 h_in: float):
        return super().exchange(mdot, P, h_in)

class Evaporator(HeatExchanger):
    def __init__(self,
                 sim: SimCycle,
                 params):

        super().__init__(sim, params.N_eva, params.UA_eva_elem, params.T_eva_air,
                    mode='eva',
                    V_elem=params.V_elem_eva)
        

    def exchange(self,
                 mdot: float,
                 P: float,
                 h_in: float):
        return super().exchange(mdot, P, h_in)
