from ..core.sim_cycle import SimCycle

class ChargeCalculator:
    def __init__(self, sim_cycle: SimCycle, target: float):
        self.sim = sim_cycle
        self.target = target

    def error(self, m_tot: float):
        return (m_tot - self.target) / self.target