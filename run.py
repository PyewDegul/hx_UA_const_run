# run.py
from hx_UA_const.cycle.cycle_model import CycleModel
from hx_UA_const.core.params import SystemParams  # 예시 경로, 실제 구조에 따라 조정
import cProfile
import pstats

params = SystemParams(
    UA_total=1000, N_cond=200, N_eva=50,
    T_cond_air=35 + 273.15, T_eva_air=27+ 273.15,
    isen_eff = 0.9, V_comp =2e-5, f_comp=50,
    DSH_target=5, DSC_target=5, tol=0.01
)
model = CycleModel(params, "REFPROP","R32")

cProfile.run('model.run()')  # Profiling the run method
results = model.run()

print(results)
