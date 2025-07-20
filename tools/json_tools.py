import importlib


class JSONAble:

    def to_dict(self):
        ret = {}
        ret["@class"] = self.__class__.__name__
        ret["@module"] = self.__class__.__module__
        return ret

    @classmethod
    def from_dict(cls, dic: dict):
        ret_m = importlib.import_module(dic.pop("@module"))
        ret_cls = getattr(ret_m, dic.pop("@class"))
        return ret_cls(**dic)
