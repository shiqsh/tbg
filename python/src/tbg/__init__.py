"""TBG 数值模型的公共接口。

导入本包只定义类和函数，不执行能带扫描，不读取或写入计算结果，也不画图。
使用示例：from tbg import BMParameters, BMModel, make_k_path, solve_bands
"""

from .bm_model import BMModel, BMParameters, make_k_path, solve_bands

__all__ = ['BMParameters', 'BMModel', 'make_k_path', 'solve_bands']
