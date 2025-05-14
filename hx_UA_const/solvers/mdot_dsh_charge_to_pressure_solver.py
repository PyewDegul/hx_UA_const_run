import scipy.optimize as opt
import numpy as np
from ..core.sim_cycle import SimCycle

from hx_UA_const.components.compressor import Compressor
from hx_UA_const.components.expansion_valve import ExpansionValve
from hx_UA_const.components.heat_exchanger import Condenser, Evaporator
from hx_UA_const.components.connector import Connector
from hx_UA_const.metrics.dsh_dsc_cal import DSHCalculator
from hx_UA_const.metrics.charge_cal import ChargeCalculator
from hx_UA_const.metrics.mdot_cal import mdotCalculator

from dataclasses import dataclass

@dataclass
class solved_cond_results:
    P_cond_sol: float
    h_comp_out: float
    s_comp_out: float
    T_comp_out: float
    h_cond_elem: np.ndarray
    s_cond_elem: np.ndarray
    T_cond_elem: np.ndarray
    h_exp_out: float
    s_exp_out: float
    T_exp_out: float
    mdot: float
    m_cond: float

@dataclass
class solved_eva_results(solved_cond_results):
    P_eva_sol: float
    h_eva_elem: np.ndarray
    s_eva_elem: np.ndarray
    T_eva_elem: np.ndarray
    m_eva: float

@dataclass
class solved_results(solved_eva_results):
    mtot: float
    DSH_sol: float

class PressureSolver_mdot_DSH_charge:
    def __init__(self,
                 sim:SimCycle,
                 params):
        self.sim = sim
        self.params = params

        self.comp = Compressor(sim, self.params)
        self.cond = Condenser(sim, self.params)
        self.exp = ExpansionValve(sim, self.params)
        self.eva = Evaporator(sim, self.params)
        self.conn = Connector(sim, self.params)
        self.dsh = DSHCalculator(sim, self.params.DSH_target)
        self.charge = ChargeCalculator(sim, self.params.charge_target)
        self.Mdot = mdotCalculator(sim)

        self.tol = params.tol
 

    def solve_cond(self, P_eva: float, T_cond_air: float):
        def cycle_mdot(P_cond):
            h_comp_out, s_comp_out, T_comp_out, mdot = self.comp.process(P_eva, P_cond)
            h_cond_elem, s_cond_elem, T_cond_elem, m_cond = self.cond.exchange(mdot, P_cond, h_comp_out)
            h_exp_out, s_exp_out, T_exp_out, _ = self.exp.process(P_eva, P_cond, h_cond_elem[-1])
            return solved_cond_results(
                P_cond_sol=P_cond,
                h_comp_out=h_comp_out,
                s_comp_out=s_comp_out,
                T_comp_out=T_comp_out,
                h_cond_elem=h_cond_elem,
                s_cond_elem=s_cond_elem,
                T_cond_elem=T_cond_elem,
                h_exp_out=h_exp_out,
                s_exp_out=s_exp_out,
                T_exp_out=T_exp_out,
                mdot=mdot,
                m_cond=m_cond
            )
        def mdot_err(P_cond):
            solved_cond_res = cycle_mdot(P_cond)
            _, _, _, mdot_exp = self.exp.process(P_eva, P_cond, solved_cond_res.h_cond_elem[-1])
            return self.Mdot.error(solved_cond_res.mdot, mdot_exp)
        
        # bisect or brentq or toms748
        # Air_Temp 부터 T_critical 사이에서 해 탐색
        # 안전 범위 0.95 * P_c
        P_cond_low = self.sim.get_single('QT_inputs', 0, T_cond_air, ('P'))
        P_cond_high = min(P_cond_low + (self.sim.P_C - P_cond_low) * 0.5, self.sim.P_C * 0.95)
        P_cond_sol = opt.brentq(mdot_err, P_cond_low, P_cond_high, xtol=self.tol)
        return cycle_mdot(P_cond_sol)
        
            
    def solve_evap(self, T_cond_air: float, T_eva_air: float):
        def cycle_DSH(P_eva):
            solved_cond = self.solve_cond(P_eva, T_cond_air)
            h_eva_elem, s_eva_elem, T_eva_elem, m_eva = self.eva.exchange(solved_cond.mdot, P_eva, solved_cond.h_exp_out)
            return solved_eva_results(
                **vars(solved_cond),
                P_eva_sol=P_eva,
                h_eva_elem=h_eva_elem,
                s_eva_elem=s_eva_elem,
                T_eva_elem=T_eva_elem,
                m_eva=m_eva
            )
        def DSH_err(P_eva):
            solved_eva = cycle_DSH(P_eva)
            return self.dsh.error(solved_eva.T_eva_elem[-1], P_eva)
        
        # bisect or brentq or toms748
        # Air_Temp 부터 T_triple 사이에서 해 탐색
        # 안전 범위 0.1 MPa
        P_eva_high = self.sim.get_single('QT_inputs', 1, T_eva_air, ('P'))
        P_eva_low = max(P_eva_high - (P_eva_high - self.sim.P_TP) * 0.5, 0.1 * 1e6)
        P_eva_sol = opt.brentq(DSH_err, P_eva_low, P_eva_high, xtol=self.tol)
        return cycle_DSH(P_eva_sol)
    
    def solve_charge(self, T_cond_air: float, T_eva_air: float):
        def cycle_charge(DSH_ass):
            self.params.DSH_target = DSH_ass
            self.comp.set_DSH(DSH_ass)
            self.dsh.set_target(DSH_ass)
            solved_eva = self.solve_evap(T_cond_air, T_eva_air)
            m_conn = self.conn.process(solved_eva.P_eva_sol, solved_eva.h_cond_elem[-1], solved_eva.h_eva_elem[-1])
            mtot = solved_eva.m_cond + solved_eva.m_eva + m_conn
            return solved_results(
                **vars(solved_eva),
                mtot=mtot,
                DSH_sol=DSH_ass
            )
        def charge_err(DSH_ass):
            solved_res = cycle_charge(DSH_ass)
            return self.charge.error(solved_res.mtot)
        
        # bisect or brentq or toms748
        DSH_low = 1e-3
        DSH_high = 20
        DSH_sol = opt.brentq(charge_err, DSH_low, DSH_high, xtol=self.tol)
        return cycle_charge(DSH_sol)

    