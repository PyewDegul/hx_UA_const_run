# cycle/cycle_model.py

from hx_UA_const.core.sim_cycle import SimCycle
from hx_UA_const.core.params import SystemParams
from hx_UA_const.components.compressor import Compressor
from hx_UA_const.components.heat_exchanger import Condenser, Evaporator
from hx_UA_const.metrics.dsh_dsc_cal import DSHCalculator, DSCCalculator
from hx_UA_const.solvers.dsh_dsc_to_pressure_solver import PressureSolver

class CycleModel:
    def __init__(self, params: SystemParams, backend: str, fluid: str):
        self.params = params
        sim = SimCycle(backend, fluid)
        self.comp = Compressor(sim, params)
        self.cond = Condenser(sim, params.N_cond,
                              params.UA_total / params.N_cond,
                              params.T_cond_air)
        self.eva = Evaporator(sim, params.N_eva,
                              params.UA_total / params.N_eva,
                              params.T_eva_air)
        self.dsh = DSHCalculator(sim, params.DSH_target)
        self.dsc = DSCCalculator(sim, params.DSC_target)
        self.solver = PressureSolver(sim, self.comp, self.cond,
                                     self.eva, self.dsh,
                                     self.dsc, params.tol)

    def run(self):
        P_cond, P_eva = self.solver.solve_cond(self.params.T_cond_air, self.params.T_eva_air)
        return {'P_cond [MPa]': P_cond / 1e6 , 'P_eva [MPa]': P_eva / 1e6}       
        