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
