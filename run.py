# run.py
from hx_UA_const.cycle.cycle_model import CycleModel, CycleModel_charge, CycleModel_mdot_charge
from hx_UA_const.core.params import SystemParams_DSH_DSC, SystemParams_DSH_charge, SystemParams_mdot_DSH_charge # 예시 경로, 실제 구조에 따라 조정

import cProfile
import pstats

'''1. DSH, DSC to Pressure Solver(2-loop)'''
params1 = SystemParams_DSH_DSC(
    UA_total=1000, N_cond=200, N_eva=50,
    T_cond_air=35 + 273.15, T_eva_air=27+ 273.15,
    isen_eff = 0.9, V_comp =2e-5, f_comp=50, 
    CA = None,
    DSH_target=5, DSC_target=5, tol=0.01
)

model1 = CycleModel(params1, "REFPROP","R32")
# cProfile.run('model1.run()')  # Profiling the run method

results = model1.run()
model1.plot()

print(results)

'''2. DSH, Charge to Pressure Solver(2-loop)'''
# DSH, Charge to Pressure Solver
params2 = SystemParams_DSH_charge(
    U_cond=1000, U_eva=1000, N_cond=200, N_eva=50,
    D_cond=8e-3, L_cond=30, D_eva=6e-3, L_eva=30, L_connect=5, 
    T_cond_air=35 + 273.15, T_eva_air=27+ 273.15,
    isen_eff = 0.7, V_comp =2e-5, f_comp=50, 
    CA = None,
    DSH_target=5, charge_target=0.32, tol=0.01
)

model2 = CycleModel_charge(params2, "REFPROP","R32")
# cProfile.run('model2.run()')  # Profiling the run method
results = model2.run()
model2.plot()
print(results)

'''3. mdot, DSH, Charge to Pressure Solver(3-loop)'''
# mdot, DSH, Charge to Pressure Solver(DSH 고정 X)
params3 = SystemParams_mdot_DSH_charge(
    U_cond=1000, U_eva=1000, N_cond=200, N_eva=50,
    D_cond=8e-3, L_cond=30, D_eva=6e-3, L_eva=30, L_connect=5, 
    T_cond_air=35 + 273.15, T_eva_air=27+ 273.15,
    isen_eff = 0.7, V_comp =2e-5, f_comp=50, 
    CA=8e-7,
    DSH_target=0.0, charge_target = 0.32, tol=0.01
)
model3 = CycleModel_mdot_charge(params3, "REFPROP","R32")
# cProfile.run('model3.run()')  # Profiling the run method
'''
results = model3.run()
model3.plot()
print(results)

'''