from ase.utils.parsemath import eval_expression
import math

param_dct = {
    "param_1": 1,
    "param_2": 4.0,
    "param_3": 478245.7586,
    "param_4": 1.58e-5,
    "param_5": -2.48757
}

expressions = [
    "3.0*param_1",
    "param_2**-2.0",
    "param_1 / param_5",
    "param_1 + param_2 / 30.0 - param_5",
    "(param_1 + 1) / param_2",
    "sqrt(param_2)",
    "fmod(param_4, param_1)",
]

solutions = [
    3*param_dct["param_1"],
    param_dct["param_2"]**-2.0,
    param_dct["param_1"] / param_dct["param_5"],
    param_dct["param_1"] + param_dct["param_2"] / 30.0 - param_dct["param_5"],
    (param_dct["param_1"] + 1) / param_dct["param_2"],
    math.sqrt(param_dct["param_2"]),
    math.fmod(param_dct["param_4"], param_dct["param_1"]),
]

for expr, soln in zip(expressions, solutions):
    assert abs(eval_expression(expr, param_dct) - soln) < 1e-13

try:
    eval_expression("99**99**99*99**99**99")
    raise RuntimeError("This should not be reached, the parser is now vulnerable to computational time based DNS attack")
except ValueError:
    pass

try:
    eval_expression("e"*10000000, dict())
    raise RuntimeError("This should not be reached, the parser is now vulnerable to memory based DNS attack")
except ValueError:
    pass

try:
    eval_expression("__import__('os').system('echo $HOME')")
    raise RuntimeError("This should not be reached, the parser can execute malicious code")
except TypeError:
    pass

