""".. Line to protect from pydocstyle D205, D400.

Metrics
=======

iCount processing statistics can be stored into instances of :py:mod:`iCount.metrics.Metrics`.

.. autoclass:: Metrics
    :members:
    :special-members: __init__

"""
import inspect


class Metrics:
    """Storge for statistics calculated during function execution."""

    def __init__(self, context=None, **kwargs):
        """
        When creating the placeholder, a process-specific context should be given.

        Parameters
        ----------
        context : str
            Context is used to indicate the process that generated the processing statistics.

        Returns
        -------
        Metrics
            Instance where processing statistics can be added or modified.

        """
        # If context is not given, determine it from calling function.
        if context is None:
            previous_frame = inspect.getouterframes(inspect.currentframe())[1]
            module = inspect.getmodulename(previous_frame[1])
            context = module + '.' + previous_frame[3]
        self.context = context
        for arg, value in kwargs.items():
            setattr(self, arg, value)

    def __repr__(self):
        """Customize repr to enable recreatin of objects with eval statement."""
        attrs = [(k, v) for k, v in self.__dict__.items() if k != 'context' and not callable(v)]
        attrs = [('context', self.context)] + attrs
        args = ', '.join('{}="{}"'.format(k, v.__repr__()) for k, v in attrs)
        return 'Metrics({:s})'.format(args)
