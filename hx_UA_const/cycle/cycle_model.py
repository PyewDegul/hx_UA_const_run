# cycle/cycle_model.py
'''1. DSH, DSC to Pressure Solver(2-loop)'''
from hx_UA_const.core.sim_cycle import SimCycle
from hx_UA_const.core.params import SystemParams_DSH_DSC
from hx_UA_const.solvers.dsh_dsc_to_pressure_solver import PressureSolver

class CycleModel:
    def __init__(self, params: SystemParams_DSH_DSC, backend: str, fluid: str):
        sim = SimCycle(backend, fluid)
        self.params = params
        self.solver = PressureSolver(sim, self.params)

    def run(self):
        res = self.solver.solve_cond(self.params.T_cond_air, self.params.T_eva_air)
        return {'P_cond [MPa]': res.P_cond_sol / 1e6 , 'P_eva [MPa]': res.P_eva_sol / 1e6}
           
'''2. DSH, Charge to Pressure Solver(2-loop)'''
from hx_UA_const.core.sim_cycle import SimCycle
from hx_UA_const.core.params import SystemParams_DSH_charge
from hx_UA_const.solvers.dsh_charge_to_pressure_solver import PressureSolver_charge

class CycleModel_charge:
    def __init__(self, params: SystemParams_DSH_charge, backend: str, fluid: str):
        sim = SimCycle(backend, fluid)
        self.params = params
        self.solver = PressureSolver_charge(sim, self.params)

    def run(self):
        res = self.solver.solve_cond(self.params.T_cond_air, self.params.T_eva_air)
        return {'P_cond [MPa]': res.P_cond_sol / 1e6 , 'P_eva [MPa]': res.P_eva_sol / 1e6} 
    
'''3. mdot, DSH, Charge to Pressure Solver(3-loop)'''
from hx_UA_const.core.sim_cycle import SimCycle
from hx_UA_const.core.params import SystemParams_mdot_DSH_charge
from hx_UA_const.solvers.mdot_dsh_charge_to_pressure_solver import PressureSolver_mdot_DSH_charge

class CycleModel_mdot_charge:
    def __init__(self, params: SystemParams_mdot_DSH_charge, backend: str, fluid: str):
        sim = SimCycle(backend, fluid)
        self.params = params
        self.solver = PressureSolver_mdot_DSH_charge(sim, self.params)

    def run(self):
        res = self.solver.solve_charge(self.params.T_cond_air, self.params.T_eva_air)
        return {'P_cond [MPa]': res.P_cond_sol / 1e6 , 'P_eva [MPa]': res.P_eva_sol / 1e6, 'DSH [K]' : res.DSH_sol}