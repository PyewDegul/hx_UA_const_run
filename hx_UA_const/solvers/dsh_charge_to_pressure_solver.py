import scipy.optimize as opt
import numpy as np
from ..core.sim_cycle import SimCycle

from hx_UA_const.components.compressor import Compressor
from hx_UA_const.components.expansion_valve import ExpansionValve
from hx_UA_const.components.heat_exchanger import Condenser, Evaporator
from hx_UA_const.components.connector import Connector
from hx_UA_const.metrics.dsh_dsc_cal import DSHCalculator
from hx_UA_const.metrics.charge_cal import ChargeCalculator

from dataclasses import dataclass

@dataclass
class solved_eva_results:
    P_eva_sol: float
    h_comp_out: float
    s_comp_out: float
    T_comp_out: float
    h_cond_elem: np.ndarray
    s_cond_elem: np.ndarray
    T_cond_elem: np.ndarray
    h_exp_out: float
    s_exp_out: float
    T_exp_out: float
    h_eva_elem: np.ndarray
    s_eva_elem: np.ndarray
    T_eva_elem: np.ndarray
    m_cond: float
    m_eva: float
    mdot: float

@dataclass
class solved_results(solved_eva_results):
    P_cond_sol: float
    mtot: float
    

class PressureSolver_charge:
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
        
        self.tol = self.params.tol


    def solve_evap(self, P_cond: float, T_eva_air: float):
        def cycle_dsh(P_eva):
            h_comp_out, s_comp_out, T_comp_out, mdot = self.comp.process(P_eva, P_cond)
            h_cond_elem, s_cond_elem, T_cond_elem, m_cond = self.cond.exchange(mdot, P_cond, h_comp_out)
            h_exp_out, s_exp_out, T_exp_out = self.exp.process(P_eva, P_cond, h_cond_elem[-1])
            h_eva_elem, s_eva_elem, T_eva_elem, m_eva = self.eva.exchange(mdot, P_eva, h_exp_out)    
            return solved_eva_results(
                P_eva_sol=P_eva,
                h_comp_out=h_comp_out,
                s_comp_out=s_comp_out,
                T_comp_out=T_comp_out,
                h_cond_elem=h_cond_elem,
                s_cond_elem=s_cond_elem,
                T_cond_elem=T_cond_elem,
                h_exp_out=h_exp_out,
                s_exp_out=s_exp_out,
                T_exp_out=T_exp_out,
                h_eva_elem=h_eva_elem,
                s_eva_elem=s_eva_elem,
                T_eva_elem=T_eva_elem,
                m_cond=m_cond,
                m_eva=m_eva,
                mdot=mdot
            )
        def dsh_err(P_eva):
            solved_eva_res = cycle_dsh(P_eva)
            return self.dsh.error(solved_eva_res.T_eva_elem[-1], P_eva)
        
        # bisect or brentq or toms748
        P_eva_high = self.sim.get_single('QT_inputs', 1, T_eva_air, ('P'))
        # P_eva_low = 0.1 * 1e6
        P_eva_low = max(P_eva_high - (P_eva_high - 0.1 * 1e6) * 0.5, 0.1 * 1e6)
        P_eva_sol = opt.brentq(dsh_err, P_eva_low, P_eva_high, xtol=self.tol)
        return cycle_dsh(P_eva_sol)
            
    def solve_cond(self, T_cond_air: float, T_eva_air: float):
        def cycle_charge(P_cond):
            solved_eva = self.solve_evap(P_cond, T_eva_air)
            m_conn = self.conn.process(solved_eva.P_eva_sol, solved_eva.h_cond_elem[-1], solved_eva.h_eva_elem[-1])
            mtot = solved_eva.m_cond + solved_eva.m_eva + m_conn
            return solved_results(
                **vars(solved_eva),  # Unpack the solved_eva dataclass
                P_cond_sol=P_cond,
                mtot=mtot
                )
        def charge_err(P_cond):
            solved_res = cycle_charge(P_cond)
            return self.charge.error(solved_res.mtot)
        
        # bisect or brentq or toms748
        P_cond_low = self.sim.get_single('QT_inputs', 0, T_cond_air, ('P'))
        # P_cond_high = self.sim.P_C
        P_cond_high = min(P_cond_low + (self.sim.P_C - P_cond_low) * 0.5, self.sim.P_C * 0.95)
        P_cond_sol = opt.brentq(charge_err, P_cond_low, P_cond_high, xtol=self.tol)
        return cycle_charge(P_cond_sol)

