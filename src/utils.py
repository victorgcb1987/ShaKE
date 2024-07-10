


def check_run(results):
    if results["returncode"] == 0:
        return "#SUCCESS: {}".format(results["command"])
    elif results["returncode"] == 99:
        return "#ALREADY_DONE: {}".format(results["command"])
    else:
        return "#FAIL: {} {}".format(results["command"], results["msg"])