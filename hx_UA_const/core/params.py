# core/params.py
import numpy as np
from dataclasses import dataclass

@dataclass
class BaseParams_DSH:
    N_cond: int
    N_eva: int
    T_cond_air: float
    T_eva_air: float
    isen_eff: float
    V_comp: float
    f_comp: float
    CA: float
    DSH_target: float
    tol: float

    def __post_init__(self):
        # Ensure that the parameters are within reasonable bounds
        assert self.N_cond > 0 and self.N_eva > 0, \
        "N_cond, N_eva must be positive"

@dataclass
class Params_Geometry_charge(BaseParams_DSH):
    # 기하 형상
    D_cond: float
    L_cond: float
    D_eva: float
    L_eva: float
    L_connect: float

    @property
    def A_elem_cond(self):
        return self.D_cond * self.L_cond * np.pi / self.N_cond
    @property
    def V_elem_cond(self):
        return self.L_cond * np.pi * self.D_cond**2 / (4*self.N_cond)
    
    @property
    def A_elem_eva(self):
        return self.D_eva * self.L_eva * np.pi / self.N_eva
    @property
    def V_elem_eva(self):
        return self.L_eva * np.pi * self.D_eva**2 / (4*self.N_eva)
        
    @property
    def V_connect(self):
        return np.pi * self.L_connect * (self.D_eva **2) / 4

    def __post_init__(self):
        for name in ("D_cond","L_cond","D_eva","L_eva","L_connect"):
            val = getattr(self, name)
            assert val > 0, f"{name} must be positive"

# dataclass for DSH, DSC
@dataclass
class SystemParams_DSH_DSC(BaseParams_DSH):
    UA_total: float
    DSC_target: float


    @property
    def UA_cond(self):
        return self.UA_total
    @property
    def UA_eva(self):
        return self.UA_total
    
    @property
    def UA_cond_elem(self):
        return self.UA_cond / self.N_cond
    @property
    def UA_eva_elem(self):
        return self.UA_eva / self.N_eva
    
    # Dummy properties for geometry
    @property
    def V_elem_cond(self):
        return 0.0
    @property
    def V_elem_eva(self):
        return 0.0
    
    def __post_init__(self):
        assert self.UA_total > 0, "UA_total must be positive"
        
# dataclass for, DSH, Charge
@dataclass
class SystemParams_DSH_charge(Params_Geometry_charge):
    U_cond: float
    U_eva: float
    charge_target: float
    
    @property
    def UA_cond_elem(self):
        return self.U_cond * self.A_elem_cond
    @property
    def UA_eva_elem(self):
        return self.U_eva * self.A_elem_eva
    
    @property
    def UA_cond(self):
        return self.UA_cond_elem * self.N_cond
    @property
    def UA_eva(self):
        return self.UA_eva_elem * self.N_eva

# dataclass for mdot, DSH, Charge
@dataclass
class SystemParams_mdot_DSH_charge(Params_Geometry_charge):
    U_cond: float
    U_eva: float
    charge_target: float
    
    @property
    def UA_cond_elem(self):
        return self.U_cond * self.A_elem_cond
    @property
    def UA_eva_elem(self):
        return self.U_eva * self.A_elem_eva
    
    @property
    def UA_cond(self):
        return self.UA_cond_elem * self.N_cond
    @property
    def UA_eva(self):
        return self.UA_eva_elem * self.N_eva

