# run.py
from hx_UA_const.cycle.cycle_model import CycleModel, CycleModel_charge
from hx_UA_const.core.params import SystemParams, SystemParams_DSH_charge # 예시 경로, 실제 구조에 따라 조정
import cProfile
import pstats

# DSH, DSC to Pressure Solver
params1 = SystemParams(
    UA_total=1000, N_cond=200, N_eva=50,
    T_cond_air=35 + 273.15, T_eva_air=27+ 273.15,
    isen_eff = 0.9, V_comp =2e-5, f_comp=50,
    DSH_target=5, DSC_target=5, tol=0.01
)

model1 = CycleModel(params1, "BICUBIC&HEOS","R32")

cProfile.run('model1.run()')  # Profiling the run method
# results = model1.run()

# print(results)



# DSH, Charge to Pressure Solver
params2 = SystemParams_DSH_charge(
    U_cond=1000, U_eva=1000, N_cond=40, N_eva=40,
    D_cond=8e-3, L_cond=30, D_eva=6e-3, L_eva=30,
    L_connect=5, T_cond_air=35 + 273.15, T_eva_air=27+ 273.15,
    isen_eff = 0.7, V_comp =2e-5, f_comp=50,
    DSH_target=5, charge_target=0.6, tol=0.01
)
model2 = CycleModel_charge(params2, "BICUBIC&HEOS","R410a")


cProfile.run('model2.run()')  # Profiling the run method
results = model2.run()

print(results)