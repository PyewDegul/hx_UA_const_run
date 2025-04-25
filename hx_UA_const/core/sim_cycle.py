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

