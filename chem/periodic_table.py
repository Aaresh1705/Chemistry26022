from sympy import Matrix, lcm


class Element:
    atomic_number = 0
    molar_mass = u = 0

    def __init__(self, n=1, grams=0, mol=0, state='', charge=0):
        self.n = n
        self.coef = 1

        self.grams = grams
        self.mol = mol

        self.state = state
        self.charge = charge

    def total_molar_mass(self):
        return self.molar_mass * self.n

    def mol(self):
        self.total_molar_mass() * self.grams

    def calculate_mol(self):
        self.mol = 1 / self.total_molar_mass() * self.grams
        return self.mol

    def calculate_gram(self):
        self.grams = self.mol * self.total_molar_mass()
        return self.grams

    def get_electrons(self):
        e = self.atomic_number * self.n + self.charge
        return e

    def get_thermal_properties(self):
        if self.state != '':
            key = (self, self.charge)
            if key in THERMO_TABLE:
                return THERMO_TABLE[key][self.state]

            key = self
            if key in THERMO_TABLE:
                return THERMO_TABLE[key][self.state]

            raise ValueError(f"No enthalpy data for {self}({self.state}) with charge {self.charge}")
        else:
            return False

    def _composition(self):
        """Return composition as a sorted tuple: [('H', 2), ('O', 1)]"""

        return tuple(sorted([(type(self).__name__, self.n)]))

    def __repr__(self):
        parts = []
        for el, count in self._composition():
            parts.append(el + (str(count) if count > 1 else ""))
        return "".join(parts)

    def __eq__(self, other):
        return isinstance(other, Element) and self._composition() == other._composition()

    def __hash__(self):
        return hash(self._composition())


class Molekyle(Element):
    def __init__(self, *elements: Element, state='', charge=0):
        self.elements = elements
        self.n = 1
        self.coef = 1

        self.grams = 0
        self.mol = 0

        self.state = state
        self.charge = charge

    def total_molar_mass(self):
        molar_mass = 0
        for element in self.elements:
            molar_mass += element.total_molar_mass()

        return molar_mass

    def _composition(self):
        """Return composition as a sorted tuple: [('H', 2), ('O', 1)]"""
        comp = {}
        for e in self.elements:
            symbol = type(e).__name__  # "H", "O", etc.
            comp[symbol] = comp.get(symbol, 0) + e.n
        return tuple(sorted(comp.items()))

    def get_electrons(self):
        e = sum([e.get_electrons() for e in self.elements]) - self.charge
        return e


class Reaction:
    def __init__(self, left: list[Element], right: list[Element]):
        self.left = left
        self.right = right

    def balance(self):
        # Collect all elements
        elements = set()
        for mol in self.left + self.right:
            for el, _ in mol._composition():
                elements.add(el)
        elements = sorted(elements)

        # Build matrix
        rows = []
        for el in elements:
            row = []
            for mol in self.left:
                comp = dict(mol._composition())
                row.append(-comp.get(el, 0))
            for mol in self.right:
                comp = dict(mol._composition())
                row.append(comp.get(el, 0))
            rows.append(row)

        M = Matrix(rows)
        nullspace = M.nullspace()
        if not nullspace:
            raise ValueError("No solution found")

        vec = nullspace[0]
        # Normalize first term to 1
        vec = vec / vec[0]
        lcm_denom = lcm([term.q for term in vec])
        coefs = [int(term * lcm_denom) for term in vec]

        for mol, coef in zip(self.left + self.right, coefs):
            mol.coef = coef

    def __repr__(self):
        def side_to_str(side):
            parts = []
            for mol in side:
                # coefficient (skip if 1)
                coef = f"{mol.coef if mol.coef != 1 else ''}"
                # state (skip if empty)
                state = f"({mol.state})" if mol.state else ""
                # charge (skip if 0)
                charge = f"[{mol.charge:+}]" if mol.charge else ""
                parts.append(f"{coef}{mol}{state}{charge}")
            return " + ".join(parts)

        return f"{side_to_str(self.left)} -> {side_to_str(self.right)}"

    def reaction_enthalpy(self):
        """Compute reaction enthalpy (ΔH) in kJ/mol, considering state and charge"""
        delta_h = 0.0
        p = 'H'

        for mol in self.right:  # products
            delta_h += mol.coef * mol.get_thermal_properties()[p]

        for mol in self.left:  # reactants
            delta_h -= mol.coef * mol.get_thermal_properties()[p]

        return delta_h

    def reaction_entropy(self):
        """Compute reaction entropy in kJ/mol"""
        delta_s = 0.0
        p = 'S'

        for mol in self.right:  # products
            delta_s += mol.coef * mol.get_thermal_properties()[p]

        for mol in self.left:  # reactants
            delta_s -= mol.coef * mol.get_thermal_properties()[p]

        return delta_s

    def reaction_energy(self):
        """Compute reaction energy in kJ/mol"""
        delta_g = 0.0
        p = 'G'

        for mol in self.right:  # products
            delta_g += mol.coef * mol.get_thermal_properties()[p]

        for mol in self.left:  # reactants
            delta_g -= mol.coef * mol.get_thermal_properties()[p]

        return delta_g

class PeriodicTable:
    class H(Element):
        atomic_number = 1
        molar_mass = u = 1.01

    class He(Element):
        atomic_number = 2
        molar_mass = u = 2.00

    class Li(Element):
        atomic_number = 3
        molar_mass = u = 6.94

    class Be(Element):
        atomic_number = 4
        molar_mass = u = 9.01

    class B(Element):
        atomic_number = 5
        molar_mass = u = 10.81

    class C(Element):
        atomic_number = 6
        molar_mass = u = 12.01

    class N(Element):
        atomic_number = 7
        molar_mass = u = 14.01

    class O(Element):
        atomic_number = 8
        molar_mass = u = 16.00

    class F(Element):
        atomic_number = 9
        molar_mass = u = 19.00

    class Ne(Element):
        atomic_number = 10
        molar_mass = u = 20.18

    class Na(Element):
        atomic_number = 11
        molar_mass = u = 22.99

    class Mg(Element):
        atomic_number = 12
        molar_mass = u = 24.31

    class Al(Element):
        atomic_number = 13
        molar_mass = u = 26.98

    class Si(Element):
        atomic_number = 14
        molar_mass = u = 28.09

    class P(Element):
        atomic_number = 15
        molar_mass = u = 30.97

    class S(Element):
        atomic_number = 16
        molar_mass = u = 32.06

    class Cl(Element):
        atomic_number = 17
        molar_mass = u = 35.45

    class Ar(Element):
        atomic_number = 18
        molar_mass = u = 39.95

    class K(Element):
        atomic_number = 19
        molar_mass = u = 39.10

    class Ca(Element):
        atomic_number = 20
        molar_mass = u = 40.08

    class Sc(Element):
        atomic_number = 21
        molar_mass = u = 44.96

    class Ti(Element):
        atomic_number = 22
        molar_mass = u = 47.87

    class V(Element):
        atomic_number = 23
        molar_mass = u = 50.94

    class Cr(Element):
        atomic_number = 24
        molar_mass = u = 52.00

    class Mn(Element):
        atomic_number = 25
        molar_mass = u = 54.94

    class Fe(Element):
        atomic_number = 26
        molar_mass = u = 55.85

    class Co(Element):
        atomic_number = 27
        molar_mass = u = 58.93

    class Ni(Element):
        atomic_number = 28
        molar_mass = u = 58.69

    class Cu(Element):
        atomic_number = 29
        molar_mass = u = 63.55

    class Zn(Element):
        atomic_number = 30
        molar_mass = u = 65.38

    class Ga(Element):
        atomic_number = 31
        molar_mass = u = 69.72

    class Ge(Element):
        atomic_number = 32
        molar_mass = u = 72.63

    class As(Element):
        atomic_number = 33
        molar_mass = u = 74.92

    class Se(Element):
        atomic_number = 34
        molar_mass = u = 78.97

    class Br(Element):
        atomic_number = 35
        molar_mass = u = 79.90

    class Kr(Element):
        atomic_number = 36
        molar_mass = u = 83.80

    class Rb(Element):
        atomic_number = 37
        molar_mass = u = 85.47

    class Sr(Element):
        atomic_number = 38
        molar_mass = u = 87.62

    class Y(Element):
        atomic_number = 39
        molar_mass = u = 88.91

    class Zr(Element):
        atomic_number = 40
        molar_mass = u = 91.22

    class Nb(Element):
        atomic_number = 41
        molar_mass = u = 92.91

    class Mo(Element):
        atomic_number = 42
        molar_mass = u = 95.95

    class Tc(Element):
        atomic_number = 43
        molar_mass = u = 98.00

    class Ru(Element):
        atomic_number = 44
        molar_mass = u = 101.07

    class Rh(Element):
        atomic_number = 45
        molar_mass = u = 102.91

    class Pb(Element):
        atomic_number = 46
        molar_mass = u = 106.42

    class Ag(Element):
        atomic_number = 47
        molar_mass = u = 107.87


H2O = Molekyle(
    PeriodicTable.H(2),
    PeriodicTable.O(1),
)
O2 = PeriodicTable.O(2)
OH = Molekyle(PeriodicTable.O(), PeriodicTable.H())
H = PeriodicTable.H()
Na = PeriodicTable.Na()
SO4 = Molekyle(PeriodicTable.O(4), PeriodicTable.S())
Na2SO4 = Molekyle(PeriodicTable.Na(2), PeriodicTable.S(), PeriodicTable.O(4))
CO2 = Molekyle(PeriodicTable.C(), PeriodicTable.O(2))
C6H12O6 = Molekyle(PeriodicTable.C(6), PeriodicTable.H(12), PeriodicTable.O(6))
SO3 = Molekyle(PeriodicTable.S(), PeriodicTable.O(3))
SO2 = Molekyle(PeriodicTable.S(), PeriodicTable.O(2))
NH3 = Molekyle(PeriodicTable.N(), PeriodicTable.H(3))
N2 = PeriodicTable.N(2)
H2 = PeriodicTable.H(2)

THERMO_TABLE = {
    O2: {
        "g": {"H": 0, "G": 0, "S": 205.2},
    },
    H2O: {
        "l": {"H": -285.8, "G": -237.1,  "S": 70.0},
        "g": {"H": -241.8, "G": -228.6,  "S": 188.8},
    },
    (OH, -1): {
        'aq': {"H": -230.0, "G": -157.2,  "S": -10.8},
    },
    (H, 1): {
        'aq': {"H": 0, "G": 0, "S": 0},
    },
    (Na, 1): {
        'aq': {"H": -240.1, "G": -261.9, "S": 59.0},
    },
    (SO4, -2): {
        'aq': {"H": -909.3, "G": -744.5, "S": 20.1},
    },
    Na2SO4: {
        's': {"H": -1_384.49, "G":  -1266.8, "S": 149.49},
    },
    CO2: {
        'g': {"H": -393.5, "G": -394.4, "S": 213.8},
    },
    C6H12O6: {
        's': {"H": -1_274.5, "G": None, "S": None},
    },
    SO3: {
        'g': {"H": -395.7, "G": -371.1, "S": 256.8},
    },
    SO2: {
        'g': {"H": -296.8, "G": -300.1, "S": 248.2},
    },
    NH3: {
        'g': {"H": -45.9, "G": -16.4, "S": 192.8},
    },
    N2: {
        'g': {"H": 0, "G": 0, "S": 191.6},
    },
    H2: {
        'g': {"H": 0, "G": 0, "S": 130.7},
    },
}

if __name__ == "__main__":
    SO3 = Molekyle(PeriodicTable.O(3), PeriodicTable.S())
    SO3.state = 'g'
    print(SO3.get_thermal_properties())
