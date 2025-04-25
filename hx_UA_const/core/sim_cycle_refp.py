
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

