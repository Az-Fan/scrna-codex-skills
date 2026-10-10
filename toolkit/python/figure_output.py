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
