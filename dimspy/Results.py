import collections

class Result(object):
    def __init__(self, result_frame, type):
        self.result_frame = result_frame
        self.type = type

class ResultsList(object):
    def __init__(self):
        self.results_dict = collections.OrderedDict()

    def append(self, result):
        self.results_dict[result.type] = result.result_frame

    def variable_limiter(self, DataAnalysisObject, limit_dict={}):
        for test, values in limit_dict.items():
                for score, l_v in values.items():
                    operator, value = l_v.split(" ")
                    df = self.results_dict[test]
                    if operator == ">":
                        columns = df[(df[score] > float(value))].index
                    elif operator == "<":
                        columns = df[(df[score] < float(value))].index
                    elif operator == "=":
                        columns = df[(df[score] == float(value))].index
                DataAnalysisObject.data_frame = DataAnalysisObject.data_frame[columns]