"""Condition interpreters."""

from .base import Interpreter, Result
from .smn1 import Smn1Interpreter
from .f8 import F8Interpreter

_REGISTRY: dict[str, Interpreter] = {
    "smn1": Smn1Interpreter(),
    "f8": F8Interpreter(),
}


def get_interpreter(interpreter_id: str) -> Interpreter | None:
    return _REGISTRY.get(interpreter_id)


__all__ = ["Interpreter", "Result", "get_interpreter"]
