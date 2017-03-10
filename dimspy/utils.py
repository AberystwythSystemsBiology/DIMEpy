def _cross_validation(data, type):
    from sklearn.model_selection import LeaveOneOut
    print type
    if type == "loo":
        loo = LeaveOneOut()
        return loo.split(data)
