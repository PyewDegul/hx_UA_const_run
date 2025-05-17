from ..core.sim_cycle import SimCycle

class mdotCalculator:
    def __init__(self, sim_cycle: SimCycle):
        self.sim = sim_cycle

    def error(self, mdot_exp: float, mdot_comp:float):
        return (mdot_comp - mdot_exp) / (mdot_exp)