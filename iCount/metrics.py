import inspect


class Metrics:

    def __init__(self, context=None, **kwargs):
        # If context is not given, determine it from calling function.
        if context is None:
            previous_frame = inspect.getouterframes(inspect.currentframe())[1]
            module = inspect.getmodulename(previous_frame[1])
            context = module + '.' + previous_frame[3]
        self.context = context
        for arg, value in kwargs.items():
            setattr(self, arg, value)

    def __repr__(self):
        return '{:s}\n{:s}'.format(self.context, self.__dict__.__repr__())
