# core/sim_cycle.py
import CoolProp.CoolProp as CP
from operator import itemgetter, attrgetter

class SimCycle:
    # Backend 설정(with Coolprop wrapper), ("REFPROP" or "HEOS", ...)
    # Fluid 설정 (ex. "R32", "R134a", "R410A", ...)
    # REFPROP 백엔드 : https://www.coolprop.org/_static/doxygen/html/class_cool_prop_1_1_r_e_f_p_r_o_p_backend.html
    # 필요하면 여기서 가져올 수 있음 : 
    def __init__(self, backend_name: str, fluid_name: str):
        self.state = CP.AbstractState(backend_name, fluid_name)
        # self.P_C = self.state.p_critical()

        # input_keys: CoolProp에서 사용하는 입력값들
        # S, H, D는 모두 질량 기준
        # input_keys : https://www.coolprop.org/_static/doxygen/html/namespace_cool_prop.html#a58e7d98861406dedb48e07f551a61efb
        self.input_keys = {
            'QT_inputs': CP.QT_INPUTS,
            'PQ_inputs': CP.PQ_INPUTS,
            'PT_inputs': CP.PT_INPUTS,
            'PS_inputs': CP.PSmass_INPUTS,
            'HQ_inputs': CP.HmassQ_INPUTS, 	
            'HP_inputs': CP.HmassP_INPUTS,
            'HT_inputs': CP.HmassT_INPUTS,   
        }
        # output_pairs: CoolProp에서 사용하는 출력값들
        # output_pairs : https://www.coolprop.org/_static/doxygen/html/namespace_cool_prop.html#a4b49eeb37210a720b188f493955d8364
        self.output_pairs = {
            'T': CP.iT,
            'P': CP.iP,
            'H': CP.iHmass,
            'S': CP.iSmass,
            'D': CP.iDmass,
            'Q': CP.iQ,
            'C': CP.iCpmass,
        }
        
    # Update the state with the given inputs
    def update(self, arg, input1, input2):
        self.state.update(self.input_keys[arg], input1, input2)

    # Get the value of a single property(float) or multiple properties(ndarray)
    def get_single_no_update(self, prop : str)-> dict:
        return self.state.keyed_output(self.output_pairs[prop])

    def get_multiple_no_update(self, props : tuple):
        vals = itemgetter(*props) (self.output_pairs)
        return tuple(map(self.state.keyed_output, vals))

    def get_single(self, arg, input1, input2, props):
        self.update(arg, input1, input2)
        return self.get_single_no_update(props)
    
    def get_multiple(self, arg, input1, input2, props):
        self.update(arg, input1, input2)
        return self.get_multiple_no_update(props)

import os
import ctypes
os.environ['RPPREFIX'] = r'C:\Program Files (x86)\REFPROP'
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

class SimCycle_Refprop:
    def __init__(self, fluid_name: str):
        # 1) 환경변수에서 REFPROP 경로 읽기, 유체 설정
        self.RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
        self.RP.SETPATHdll(os.environ['RPPREFIX'])
        self.RP.SETFLUIDSdll(fluid_name)
        self.z = (ctypes.c_double * 20)(1.0, *([0.0]*19)) # composition 순수 물질
        self.iFlag = 1 # SI 질량 단위계 사용

        # 2) 인풋(ab 코드) 매핑
        self.input_map = {
            'QT_inputs': 'QT', 'PQ_inputs': 'PQ',
            'PT_inputs': 'TP', 'PS_inputs': 'PS',
            'HQ_inputs': 'HQ', 'HP_inputs': 'HP',
            'HT_inputs': 'HT',
        }
        # 3) Output 배열 인덱스(대소문자 주의)
        self.prop_index = {
            'T':  0,  'P':  1, 'D':  2, 'Dl': 3, 'Dv': 4,
            'q':  5,  'e':  6, 'h':  7, 's':  8, 'Cv': 9,
            'C': 10, 'w': 11
        }

    def update(self, arg: str, a: float, b: float, props: tuple[str, ...]):
        self.out = self.RP.ABFLSHdll(ab = self.input_map[arg], a=a, b=b, z=self.z, iFlag=self.iFlag)
        return attrgetter(*props) (self.out)


# SimCycle, SimCycle_Refprop의 비교. 10배 정도 CoolProp이 빠름

'''
import time

cycle_py = SimCycle('REFPROP', 'R32')
cycle_rp = SimCycle_Refprop('R32')

arg= 'PQ_inputs'
a, b  = 101.325, 0.0
a_prop = a * 1000
props = ('T')


# 4) 벤치마크 반복 횟수
N = 100_000

    # 5) CoolProp HEOS 측정
t0 = time.perf_counter()
for _ in range(N):
    cycle_py.update(arg, a, b, props)
t_py = time.perf_counter() - t0

    # 6) REFPROP 측정
t0 = time.perf_counter()
for _ in range(N):
    cycle_rp.update(arg, a, b, props)
t_rp = time.perf_counter() - t0

    # 7) 결과 출력
print(f"CoolProp (REFPROP)  : {t_py:.3f} s for {N} cycles")
print(f"REFPROP (ABFLSHdll): {t_rp:.3f} s for {N} cycles")

'''