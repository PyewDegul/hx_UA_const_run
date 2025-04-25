# cycle/cycle_model.py
from hx_UA_const.core.sim_cycle import SimCycle
from hx_UA_const.core.params import SystemParams

from hx_UA_const.components.compressor import Compressor
from hx_UA_const.components.heat_exchanger import Condenser, Evaporator, Condenser_charge, Evaporator_charge
from hx_UA_const.metrics.dsh_dsc_cal import DSHCalculator, DSCCalculator
from hx_UA_const.solvers.dsh_dsc_to_pressure_solver import PressureSolver


class CycleModel:
    def __init__(self, params: SystemParams, backend: str, fluid: str):
        self.params = params
        sim = SimCycle(backend, fluid)
        self.comp = Compressor(sim, params)
        self.cond = Condenser(sim, params.N_cond,
                              params.UA_cond / params.N_cond,
                              params.T_cond_air)
        self.eva = Evaporator(sim, params.N_eva,
                              params.UA_eva / params.N_eva,
                              params.T_eva_air)
        self.dsh = DSHCalculator(sim, params.DSH_target)
        self.dsc = DSCCalculator(sim, params.DSC_target)
        self.solver = PressureSolver(sim, self.comp, self.cond,
                                     self.eva, self.dsh,
                                     self.dsc, params.tol)

    def run(self):
        P_cond, P_eva = self.solver.solve_cond(self.params.T_cond_air, self.params.T_eva_air)
        return {'P_cond [MPa]': P_cond / 1e6 , 'P_eva [MPa]': P_eva / 1e6}       
    
from hx_UA_const.core.sim_cycle import SimCycle
from hx_UA_const.core.params import SystemParams_DSH_charge

from hx_UA_const.components.compressor import Compressor
from hx_UA_const.components.heat_exchanger import Condenser_charge, Evaporator_charge
from hx_UA_const.components.connector import Connector
from hx_UA_const.metrics.dsh_dsc_cal import DSHCalculator
from hx_UA_const.metrics.charge_cal import ChargeCalculator
from hx_UA_const.solvers.dsh_charge_to_pressure_solver import PressureSolver_charge

class CycleModel_charge:
    def __init__(self, params: SystemParams_DSH_charge, backend: str, fluid: str):
        self.params = params
        sim = SimCycle(backend, fluid)
        self.comp = Compressor(sim, params)
        self.cond = Condenser_charge(sim, params.N_cond,
                              params.UA_cond / params.N_cond,
                              params.T_cond_air, params.V_elem_cond)
        self.eva = Evaporator_charge(sim, params.N_eva,
                              params.UA_eva / params.N_eva,
                              params.T_eva_air, params.V_elem_eva)
        self.conn = Connector(sim, params)
        self.dsh = DSHCalculator(sim, params.DSH_target)
        self.charge = ChargeCalculator(sim, params.charge_target)
        self.solver = PressureSolver_charge(sim, self.comp, self.cond,
                                     self.eva, self.conn, self.dsh,
                                     self.charge, params.tol)

    def run(self):
        P_cond, P_eva = self.solver.solve_cond(self.params.T_cond_air, self.params.T_eva_air)
        return {'P_cond [MPa]': P_cond / 1e6 , 'P_eva [MPa]': P_eva / 1e6} 