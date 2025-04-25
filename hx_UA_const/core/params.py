# core/params.py
import numpy as np
from dataclasses import dataclass
# dataclass for DSH, DSC
@dataclass
class SystemParams:
    UA_total: float
    N_cond: int
    N_eva: int
    T_cond_air: float
    T_eva_air: float
    isen_eff: float
    V_comp: float
    f_comp: float
    DSH_target: float
    DSC_target: float
    tol: float

    def __post_init__(self):
        # Ensure that the parameters are within reasonable bounds
        assert self.UA_total > 0 and self.N_cond > 0 and self.N_eva > 0, \
        "UA_total, N_cond, N_eva must be positive"
        self.UA_cond = self.UA_total
        self.UA_eva =  self.UA_total

# dataclass for, DSH, Charge
@dataclass
class SystemParams_DSH_charge:
    U_cond: float
    U_eva: float
    N_cond: int
    N_eva: int
    # 응축기 기하 형상
    D_cond: float
    L_cond: float
    
    # 증발기 기하 형상
    D_eva: float
    L_eva: float
    
    # UA_tot 구하기
    L_connect: float
    T_cond_air: float
    T_eva_air: float
    isen_eff: float
    V_comp: float
    f_comp: float
    DSH_target: float
    charge_target: float
    tol: float

    def __post_init__(self):
        # 응축기 기하 형상
        self.A_elem_cond = self.D_cond * self.L_cond * np.pi / self.N_cond
        self.V_elem_cond = self.L_cond * np.pi * self.D_cond**2 / (4*self.N_cond) 
        # 증발기 기하 형상
        self.A_elem_eva = self.D_eva * self.L_eva * np.pi / self.N_eva
        self.V_elem_eva = self.L_eva * np.pi * self.D_eva**2 / (4*self.N_eva)
        # UA 값 계산
        self.UA_cond = self.U_cond * self.A_elem_cond * self.N_cond
        self.UA_eva = self.U_eva * self.A_elem_eva * self.N_eva
        # connector 부피
        self.V_connect = np.pi * self.L_connect * (self.D_eva **2) / 4
