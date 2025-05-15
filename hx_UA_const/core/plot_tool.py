import matplotlib
matplotlib.use('tkagg')

from matplotlib import font_manager
from matplotlib import ticker
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt

class Graph_class():
    def __init__(self):
        self.font_name = font_manager.FontProperties(fname = "c:/Windows/Fonts/Arial.ttf").get_name()
    
    def setPlotStyle(self, ax = [], ax2 = [], title = 'title', x_label = 'x_label', x_dp = 1, x_tick = 0, y_label = 'y_label', y_dp = 1, y_tick = 0, grid = 'off'):
        if grid == 'off':
            plt.rcParams['axes.grid'] = 0
        else:
            plt.rcParams['axes.grid'] = 1

        font_size = 4                               # Font 크기
        line_width = 0.65                           # 선 두께    

        ax.tick_params(width=line_width, labelsize=4.5, length=3)

        plt.rcParams['font.family'] = self.font_name    # Font 설정
        plt.rcParams['font.size'] = font_size           # Font 크기
        plt.rcParams['legend.fontsize'] = font_size     # 범례 Font 크기
        plt.rcParams['axes.xmargin'] = 0                # x축 여백
        plt.rcParams['axes.ymargin'] = 0.1                # y축 여백
        plt.rcParams['axes.formatter.useoffset'] = 'True'   # 오프셋 사용 여부
        plt.rcParams['mathtext.fontset'] = 'dejavuserif'    # 수식 Font 설정
        plt.rcParams['figure.dpi'] = 500                    # 그래프 해상도
        plt.rcParams['savefig.dpi'] = 500                   # 그래프 저장 시 해상도  
        plt.rcParams['xtick.major.pad'] = 2                 # x축 눈금과 숫자 사이 간격 
        plt.rcParams['ytick.major.pad'] = 2                 # y축 눈금과 숫자 사이 간격  
        plt.rcParams['xtick.direction'] = 'in'              # x축 눈금 방향
        plt.rcParams['ytick.direction'] = 'in'              # y축 눈금 방향
        plt.rcParams['axes.titlepad'] = 5                   # 제목과 그래프 사이 간격
        plt.rcParams['axes.linewidth'] = line_width         # 축 두께
        # plt.rcParams.keys()                               # 모든 설정 확인

        plt.title(title, weight = 'bold', size = 6)         # 그래프 제목 설정
