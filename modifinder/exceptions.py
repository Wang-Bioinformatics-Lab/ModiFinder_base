"""
**********
Exceptions
**********

Base exceptions and errors for ModiFinder.
"""

__all__ = [
    "ModiFinderException",
    "ModiFinderError",
]

class ModiFinderException(Exception):
    """Base class for exceptions in ModiFinder."""


class ModiFinderError(ModiFinderException):
    """Exception for a serious error in NetworkX"""

class ModiFinderNetworkError(ModiFinderError):
    """Exception for a serious error in NetworkX"""
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

class ModiFinderNotImplementedError(ModiFinderError):
    """Exception for a serious error in NetworkX"""
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class ModiFinderNotSolvableError(ModiFinderError):
    """Exception for a serious error in NetworkX"""
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)