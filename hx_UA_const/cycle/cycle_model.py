# cycle/cycle_model.py
from hx_UA_const.cycle.plot_model import PlotModel

'''1. DSH, DSC to Pressure Solver(2-loop)'''
from hx_UA_const.core.sim_cycle import SimCycle
from hx_UA_const.core.params import SystemParams_DSH_DSC
from hx_UA_const.solvers.dsh_dsc_to_pressure_solver import PressureSolver

from dataclasses import dataclass, asdict
from typing import Annotated, get_type_hints

@dataclass
class Cycle_model_res:
    P_cond : Annotated[float, 'MPa']
    P_eva : Annotated[float, 'MPa']
    Q_cond : Annotated[float, 'kW']
    Q_eva : Annotated[float, 'kW']
    W_comp : Annotated[float, 'kW']
    COP_H : float
    COP_R : float
    mtot : Annotated[float, 'kg']
    mdot : Annotated[float, 'kg/s']
    DSH : Annotated[float, 'K']
    DSC : Annotated[float, 'K']

    def to_dict(self, key_with_unit: bool = False)->dict:
        raw = asdict(self)
        if not key_with_unit:
            return raw
        
        hints = get_type_hints(self.__class__, include_extras=True)
        out = {}
        for key, val in raw.items():
            ann = hints.get(key)
            if getattr(ann, '__metadata__', None):
                unit = ann.__metadata__[0]
                out[f"{key} [{unit}]"] = val
            else:
                out[key] = val
        return out
    
'''1. DSH, DSC to Pressure Solver(2-loop)'''

class CycleModel:
    def __init__(self, params: SystemParams_DSH_DSC, backend: str, fluid: str):
        self.sim = SimCycle(backend, fluid)
        self.backend = backend
        self.fluid = fluid

        self.params = params
        self.solver = PressureSolver(self.sim, self.params)
        self.res = self.solver.solve_cond(self.params.T_cond_air, self.params.T_eva_air)
    
    def calculate(self):
        P_cond = self.res.P_cond_sol
        P_eva = self.res.P_eva_sol
        Q_cond =  self.res.mdot * (self.res.h_cond_elem[0] - self.res.h_cond_elem[-1])
        Q_eva = self.res.mdot * (self.res.h_eva_elem[-1] - self.res.h_eva_elem[0])
        T_eva_vap = self.sim.get_single('PQ_inputs', self.res.P_eva_sol, 1, 'T')
        # DSH_target 고정
        W_comp = self.res.mdot * (self.res.h_comp_out - self.sim.get_single('PT_inputs', self.res.P_eva_sol, T_eva_vap + self.params.DSH_target, 'H'))
        COP_H = Q_cond / W_comp
        COP_R = Q_eva / W_comp
        mtot = -999999
        mdot = self.res.mdot
        DSH = self.params.DSH_target
        DSC = self.params.DSC_target
        return Cycle_model_res(
            P_cond=P_cond / 1e6,
            P_eva=P_eva / 1e6,
            Q_cond=Q_cond / 1e3,
            Q_eva=Q_eva / 1e3,
            W_comp=W_comp / 1e3,
            COP_H=COP_H,
            COP_R=COP_R,
            mtot=mtot,
            mdot=mdot,
            DSH=DSH,
            DSC=DSC
        )

    def plot(self):
        plot_cycle = PlotModel(self.sim, self.params, self.res)
        plot_cycle.plot_TS()
        plot_cycle.plot_PH()


'''2. DSH, Charge to Pressure Solver(2-loop)'''
from hx_UA_const.core.sim_cycle import SimCycle
from hx_UA_const.core.params import SystemParams_DSH_charge
from hx_UA_const.solvers.dsh_charge_to_pressure_solver import PressureSolver_charge

class CycleModel_charge:
    def __init__(self, params: SystemParams_DSH_charge, backend: str, fluid: str):
        self.sim = SimCycle(backend, fluid)
        self.backend = backend
        self.fluid = fluid

        self.params = params
        self.solver = PressureSolver_charge(self.sim, self.params)

        self.res = self.solver.solve_cond(self.params.T_cond_air, self.params.T_eva_air)
    
    def calculate(self):
        P_cond = self.res.P_cond_sol
        P_eva = self.res.P_eva_sol
        Q_cond =  self.res.mdot * (self.res.h_cond_elem[0] - self.res.h_cond_elem[-1])
        Q_eva = self.res.mdot * (self.res.h_eva_elem[-1] - self.res.h_eva_elem[0])
        T_eva_vap = self.sim.get_single('PQ_inputs', self.res.P_eva_sol, 1, 'T')
        # DSH_target 고정
        W_comp = self.res.mdot * (self.res.h_comp_out - self.sim.get_single('PT_inputs', self.res.P_eva_sol, T_eva_vap + self.params.DSH_target, 'H'))
        COP_H = Q_cond / W_comp
        COP_R = Q_eva / W_comp
        mtot = self.res.mtot
        mdot = self.res.mdot
        DSH = self.params.DSH_target
        # DSC는 직접 계산 必要 (Q>0의 경우)
        Q_cond_liq = self.sim.get_single('HP_inputs', self.res.h_cond_elem[-1], P_cond, 'Q')
        if Q_cond_liq > 0:
            DSC = -999999
        else:
            T_cond_liq = self.sim.get_single('PQ_inputs', self.res.P_cond_sol, 0, 'T')
            DSC = abs(self.res.T_cond_elem[-1] - T_cond_liq)
        return Cycle_model_res(
            P_cond=P_cond / 1e6,
            P_eva=P_eva / 1e6,
            Q_cond=Q_cond / 1e3,
            Q_eva=Q_eva / 1e3,
            W_comp=W_comp / 1e3,
            COP_H=COP_H,
            COP_R=COP_R,
            mtot=mtot,
            mdot=mdot,
            DSH=DSH,
            DSC=DSC
        )
    def plot(self):
        plot_cycle = PlotModel(self.sim, self.params, self.res)
        plot_cycle.plot_TS()
        plot_cycle.plot_PH()
        
'''3. mdot, DSH, Charge to Pressure Solver(3-loop)'''
from hx_UA_const.core.sim_cycle import SimCycle
from hx_UA_const.core.params import SystemParams_mdot_DSH_charge
from hx_UA_const.solvers.mdot_dsh_charge_to_pressure_solver import PressureSolver_mdot_DSH_charge

class CycleModel_mdot_charge:
    def __init__(self, params: SystemParams_mdot_DSH_charge, backend: str, fluid: str):
        self.sim = SimCycle(backend, fluid)
        self.backend = backend
        self.fluid = fluid
        
        self.params = params
        self.solver = PressureSolver_mdot_DSH_charge(self.sim, self.params)

        self.res = self.solver.solve_charge(self.params.T_cond_air, self.params.T_eva_air)

    
    def calculate(self):
        P_cond = self.res.P_cond_sol
        P_eva = self.res.P_eva_sol
        Q_cond =  self.res.mdot * (self.res.h_cond_elem[0] - self.res.h_cond_elem[-1])
        Q_eva = self.res.mdot * (self.res.h_eva_elem[-1] - self.res.h_eva_elem[0])
        T_eva_vap = self.sim.get_single('PQ_inputs', self.res.P_eva_sol, 1, 'T')
        # DSH_target 변동 -> self.res의 DSH_sol로 계산 必要
        W_comp = self.res.mdot * (self.res.h_comp_out - self.sim.get_single('PT_inputs', self.res.P_eva_sol, T_eva_vap + self.res.DSH_sol, 'H'))
        COP_H = Q_cond / W_comp
        COP_R = Q_eva / W_comp
        mtot = self.res.mtot
        mdot = self.res.mdot
        # DSH_target 변동 -> self.res의 DSH_sol로 계산 必要
        DSH = self.res.DSH_sol
        # DSC는 직접 계산 必要 (Q>0의 경우)
        Q_cond_liq = self.sim.get_single('HP_inputs', self.res.h_cond_elem[-1], P_cond, 'Q')
        if Q_cond_liq > 0:
            DSC = -999999
        else:
            T_cond_liq = self.sim.get_single('PQ_inputs', self.res.P_cond_sol, 0, 'T')
            DSC = abs(self.res.T_cond_elem[-1] - T_cond_liq)     
        return Cycle_model_res(
            P_cond=P_cond / 1e6,
            P_eva=P_eva / 1e6,
            Q_cond=Q_cond / 1e3,
            Q_eva=Q_eva / 1e3,
            W_comp=W_comp / 1e3,
            COP_H=COP_H,
            COP_R=COP_R,
            mtot=mtot,
            mdot=mdot,
            DSH=DSH,
            DSC=DSC
        )

    def plot(self):
        plot_cycle = PlotModel(self.sim, self.params, self.res)
        plot_cycle.plot_TS()
        plot_cycle.plot_PH()
