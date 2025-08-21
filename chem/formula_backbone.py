import sympy as sp


class Variable:
    def __init__(self, symbol: sp.Symbol, value: float | int | None):
        self.symbol = symbol
        self.value = None
        self.status = False

        self.update_value(value)

    def update_value(self, new_value):
        self.value = new_value
        if new_value is not None:
            self.status = True

    def __mul__(self, other):
        if isinstance(other, Variable):
            return self.symbol * other.symbol
        return self.symbol * other

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, Variable):
            return self.symbol / other.symbol
        return self.symbol / other

    def __rtruediv__(self, other):
        if isinstance(other, Variable):
            return other.symbol / self.symbol
        return other / self.symbol

    def __sub__(self, other):
        if isinstance(other, Variable):
            return self.symbol - other.symbol
        return self.symbol - other

    def __rsub__(self, other):
        if isinstance(other, Variable):
            return other.symbol - self.symbol
        return other - self.symbol

    def __add__(self, other):
        if isinstance(other, Variable):
            return self.symbol + other.symbol
        return self.symbol + other

    def __radd__(self, other):
        if isinstance(other, Variable):
            return other.symbol + self.symbol
        return other + self.symbol

    def __pow__(self, power):
        if isinstance(power, Variable):
            return self.symbol ** power.symbol
        return self.symbol ** power

    def __rpow__(self, base):
        if isinstance(base, Variable):
            return base.symbol ** self.symbol
        return base ** self.symbol

    def __repr__(self):
        return f"Variable('{self.symbol}', {self.value})"


class Formula:
    formula: sp.Symbol
    _vars: dict[str, Variable]

    def __init__(self, *args):
        names = list(self.__annotations__.keys())

        symbols = sp.symbols(' '.join(names))
        self._vars = {name: Variable(sym, val) for name, sym, val in
                      zip(names, symbols, args)}

        self.symbols = symbols

    def __getattr__(self, name):
        if name in self._vars:
            return self._vars[name]
        raise AttributeError(f"This formula has no attribute '{name}'")

    def __setattr__(self, name, value):
        if name.startswith("_") or name not in getattr(self, "_vars", {}):
            super().__setattr__(name, value)
        else:
            self._vars[name].update_value(value)

    def solve(self, variable: Variable):
        var = variable.symbol
        formula = sp.solve(self.formula, var)
        return formula[0]

    def compute(self):
        # Identify which variables are missing
        missing_vars = [var for var in self._vars.values() if not var.status]
        if len(missing_vars) > 1:
            raise ValueError("Cannot compute: more than one variable is missing.")
        elif len(missing_vars) == 1:
            # Solve for the missing variable
            target = missing_vars[0]
            solved_expr = self.solve(target)
            subs_dict = {v.symbol: v.value for v in self._vars.values() if v.status}
            result = solved_expr.subs(subs_dict).simplify()
            return result.evalf()

        # If no variables are missing, just substitute into the given formula
        subs_dict = {v.symbol: v.value for v in self._vars.values()}
        result = self.formula.subs(subs_dict).simplify()
        return result.evalf()
