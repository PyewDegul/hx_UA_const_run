# core/params.py

from dataclasses import dataclass

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
