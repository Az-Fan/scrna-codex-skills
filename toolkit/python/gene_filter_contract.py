"""校验基因检出门槛，不依赖分析环境。"""
import math


def validate_thresholds(min_cells=1, min_fraction=0, min_cells_field="gene_filter.min_cells"):
    if type(min_cells) is not int or min_cells < 1:
        raise ValueError(min_cells_field + ' must be an integer >= 1')
    if isinstance(min_fraction, bool) or not isinstance(min_fraction, (int, float)) or not math.isfinite(min_fraction) or not 0 <= min_fraction <= 1:
        raise ValueError('gene_filter.min_fraction must be finite and in [0, 1]')

