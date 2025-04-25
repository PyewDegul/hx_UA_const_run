# components/heat_exchanger.py
import numpy as np
from ..core.sim_cycle import SimCycle

class HeatExchanger:
    def __init__(self, sim : SimCycle, N:int, UA_elem: float, T_air:float):
        self.sim = sim
        self.N = N
        self.UA = UA_elem
        self.T_air = T_air

    def _elementwise(self, mdot: float, P: float, h_in: float, mode: str):
        h_elem = np.empty(self.N + 1)
        T_elem = np.empty(self.N)
        q_elem = np.empty(self.N)
        Cp_elem = np.empty(self.N)
        h_elem[0] = h_in

        # iteration
        for i in range(self.N):
            # Element에서 온도, 건도, 비열 계산
            T_elem[i], q_elem[i] = self.sim.get_multiple('HP_inputs', h_elem[i], P, ('T', 'Q'))
            
            # 건도에 따른 열전달량 계산
            T_delta = abs(T_elem[i] - self.T_air)
            if 0 < q_elem[i] < 1:
                Q = self.UA * T_delta
            else:
                Cp_elem[i] = self.sim.get_single_no_update('C')
                eps = 1 - np.exp(- self.UA / (mdot * Cp_elem[i]))
                Q = eps * mdot * Cp_elem[i] * T_delta
            
            # Cond, Evap에서 엔탈피 계산
            if mode == 'cond':
                h_elem[i+1] = h_elem[i] - Q / mdot
            else:
                h_elem[i+1] = h_elem[i] + Q / mdot
        
        
        # Outlet 엔탈피, 온도 계산

        h_out = h_elem[-1]
        T_out = self.sim.get_single('HP_inputs', h_out, P, ('T'))

        return h_out, T_out


class HeatExchanger_charge:
    def __init__(self, sim : SimCycle, N:int, UA_elem: float, T_air:float, V_elem: float):
        self.sim = sim
        self.N = N
        self.UA = UA_elem
        self.V_elem = V_elem
        self.T_air = T_air
        

    def _elementwise(self, mdot: float, P: float, h_in: float, mode: str):
        h_elem = np.empty(self.N + 1)
        T_elem = np.empty(self.N)
        q_elem = np.empty(self.N)
        D_elem = np.empty(self.N)
        m_elem = np.empty(self.N)
        Cp_elem = np.empty(self.N)
        h_elem[0] = h_in

        # iteration
        for i in range(self.N):
            # Element에서 온도, 건도, 비열 계산
            T_elem[i], q_elem[i], D_elem[i] = self.sim.get_multiple('HP_inputs', h_elem[i], P, ('T', 'Q', 'D'))

            # 건도에 따른 열전달량 계산
            T_delta = abs(T_elem[i] - self.T_air)
            if 0 < q_elem[i] < 1:
                Q = self.UA * T_delta
            else:
                Cp_elem[i] = self.sim.get_single_no_update('C')
                eps = 1 - np.exp(- self.UA / (mdot * Cp_elem[i]))
                Q = eps * mdot * Cp_elem[i] * T_delta

            # Cond, Evap에서 엔탈피 계산
            if mode == 'cond':
                h_elem[i+1] = h_elem[i] - Q / mdot
            else:
                h_elem[i+1] = h_elem[i] + Q / mdot

        # Outlet 엔탈피, 온도, 내부 charge 계산
        m_elem = self.V_elem * D_elem
        m_tot = np.sum(m_elem)
        h_out = h_elem[-1]
        T_out = self.sim.get_single('HP_inputs', h_out, P, ('T'))
      
        return h_out, T_out, m_tot



class Condenser(HeatExchanger):
    def exchange(self, mdot, P, h_in):
        return self._elementwise(mdot, P, h_in, 'cond')
    
class Evaporator(HeatExchanger):
    def exchange(self, mdot, P, h_in):
        return self._elementwise(mdot, P, h_in, 'eva')
    

class Condenser_charge(HeatExchanger_charge):
    def exchange(self, mdot, P, h_in):
        return self._elementwise(mdot, P, h_in, 'cond')
    
class Evaporator_charge(HeatExchanger_charge):
    def exchange(self, mdot, P, h_in):
        return self._elementwise(mdot, P, h_in, 'eva')
    
    
