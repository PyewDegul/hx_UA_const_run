# metrics/dsh_dsc_cal.py
from ..core.sim_cycle import SimCycle

class DSHCalculator:
    def __init__(self, sim_cycle: SimCycle, target: float):
        self.sim = sim_cycle
        self.target = target

    def error(self, T_eva_out: float, P_eva: float):
        T_sat = self.sim.get_single('PQ_inputs', P_eva, 1, ('T'))
        DSH_cal = T_eva_out - T_sat
        
        return (DSH_cal - self.target) / self.target    

class DSCCalculator:
    def __init__(self, sim_cycle: SimCycle, target: float):
        self.sim = sim_cycle
        self.target = target

    def error(self, T_cond_out: float, P_cond: float):
        T_sat = self.sim.get_single('PQ_inputs', P_cond, 0, ('T'))
        DSC_cal = T_sat - T_cond_out
        return (self.target - DSC_cal) / self.target

class ChargeCalculator:
    def __init__(self, sim_cycle: SimCycle, target: float):
        self.sim = sim_cycle
        self.target = target

    def error(self, m_tot: float):
        return (m_tot - self.target) / self.target