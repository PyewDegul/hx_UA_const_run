import numpy as np

class Connector:
    def __init__(self, sim, params):
        self.sim = sim
        self.p = params
        self.V_tot = self.p.V_connect

    def process(self, P_eva: float, h_eva_in:float, h_eva_out:float) -> tuple[float, float]:
        # 증발기 전후 밀도 구하기
        rho_eva_rear = self.sim.get_single('HP_inputs', h_eva_in, P_eva, ('D'))
        rho_eva_front = self.sim.get_single('HP_inputs', h_eva_out, P_eva, ('D'))
        # Connector의 유량 계산
        m_tot = rho_eva_rear * self.V_tot + rho_eva_front * self.V_tot
        
        return m_tot