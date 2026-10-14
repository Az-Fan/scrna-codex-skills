"""A single scientific figure format per run, shared by all plotting adapters."""
def figure_format(config):
    def get(path):
        value = config
        for key in path.split('.'):
            if not isinstance(value, dict):
                return None
            value = value.get(key)
        return value
    choices = []
    for field in ('output.figure_format', 'plots.figure_format', 'enrichment.plot_format', 'visualization.figure_format'):
        value = get(field)
        if value is None:
            continue
        if not isinstance(value, str) or value.lower() not in {'pdf', 'png'}:
            raise ValueError(field + ' must select exactly one format: pdf or png')
        choices.append(value.lower())
    preview = get('output.preview_png')
    if preview is not None:
        if not isinstance(preview, bool):
            raise ValueError('output.preview_png must be true or false')
        if preview:
            choices.append('png')
    if len(set(choices)) > 1:
        raise ValueError('Conflicting figure formats; select one format for the entire run')
    return choices[0] if choices else 'pdf'


def umap_marker_area(n_cells, override=None):
    """Matplotlib/Scanpy scatter area in points squared; not R linear size."""
    import math
    if override is not None:
        if isinstance(override, bool) or not isinstance(override, (int, float)) or not math.isfinite(override) or override <= 0:
            raise ValueError('UMAP marker area must be a positive finite number')
        return override
    if n_cells < 1:
        raise ValueError('UMAP requires at least one cell')
    return min(36., max(.5, 400. / math.sqrt(n_cells)))


def categorical_palette(labels, configured=None):
    import colorsys
    import hashlib
    import warnings
    from matplotlib.colors import to_rgba
    labels = sorted(set(str(x) for x in labels))
    configured = configured or {}
    if not isinstance(configured, dict):
        raise ValueError('Categorical colours must be a label-to-colour mapping')
    signatures = [to_rgba(value) for value in configured.values()]
    if len(signatures) != len(set(signatures)):
        raise ValueError('Duplicate categorical colours, including equivalent colour aliases')
    if len(labels) > 20:
        warnings.warn('More than 20 categories: inspect colour distinguishability and retain labels', UserWarning)
    result = {key: value for key, value in configured.items() if key in labels}
    used = {to_rgba(value) for value in result.values()}
    for label in labels:
        if label in result:
            continue
        hue = int(hashlib.sha256(label.encode()).hexdigest()[:12], 16) / float(16**12)
        for attempt in range(10001):
            rgb = colorsys.hsv_to_rgb((hue + attempt * .1031) % 1., .65, .82 - .03 * (attempt % 3))
            color = '#' + ''.join(f'{round(x * 255):02X}' for x in rgb)
            if to_rgba(color) not in used:
                break
        else:
            raise ValueError('Unable to construct unique categorical colours')
        result[label] = color
        used.add(to_rgba(color))
    return result
