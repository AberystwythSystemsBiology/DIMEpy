def _cross_validation(data, type):
    from sklearn.model_selection import LeaveOneOut
    if type == "loo":
        loo = LeaveOneOut()
        return loo.split(data)
