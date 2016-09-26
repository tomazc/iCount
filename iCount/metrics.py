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
        attrs = [(k, v) for k, v in self.__dict__.items() if k != 'context' and not callable(v)]
        attrs = [('context', self.context)] + attrs
        args = ', '.join('{}="{}"'.format(k, v.__repr__()) for k, v in attrs)
        return 'Metrics({:s})'.format(args)
