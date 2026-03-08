def cli_flag(params_dict, key, flag_name, is_bool=False):
    """Return CLI flag list for a nullable config param, or [] if null.

    For non-bool params: returns [flag_name, str(value)] if value is not None.
    For bool params: returns [flag_name] if value is truthy, [] otherwise.
    """
    val = params_dict.get(key)
    if val is None:
        return []
    if is_bool:
        return [flag_name] if val else []
    return [flag_name, str(val)]


def cli_flags_str(params_dict, *specs):
    """Build a shell-ready string of optional CLI flags.

    Each spec is (key, flag_name) or (key, flag_name, True) for bool flags.
    """
    parts = []
    for spec in specs:
        is_bool = len(spec) > 2 and spec[2]
        parts.extend(cli_flag(params_dict, spec[0], spec[1], is_bool=is_bool))
    return " ".join(parts)
