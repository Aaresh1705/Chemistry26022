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
    class H:
        atomic_number = 1
        symbol = "H"
        name = "Hydrogen"

        category = "nonmetal"
        electron_configuration = "1s1"

        group = 1
        period = 1

        molar_mass = 1.01
        electronegativity = 2.20
        oxidation_states = [-1, +1]

    class He:
        atomic_number = 2
        symbol = "He"
        name = "Helium"

        category = "noble gas"
        electron_configuration = "1s2"

        group = 18
        period = 1

        molar_mass = 4.00
        electronegativity = None
        oxidation_states = [0]

    class Li:
        atomic_number = 3
        symbol = "Li"
        name = "Lithium"

        category = "alkali metal"
        electron_configuration = "[He] 2s1"

        group = 1
        period = 2

        molar_mass = 6.94
        electronegativity = 0.98
        oxidation_states = [+1]

    class Be:
        atomic_number = 4
        symbol = "Be"
        name = "Beryllium"

        category = "alkaline earth metal"
        electron_configuration = "[He] 2s2"

        group = 2
        period = 2

        molar_mass = 9.01
        electronegativity = 1.57
        oxidation_states = [+2]

    class B:
        atomic_number = 5
        symbol = "B"
        name = "Boron"

        category = "metalloid"
        electron_configuration = "[He] 2s2 2p1"

        group = 13
        period = 2

        molar_mass = 10.81
        electronegativity = 2.04
        oxidation_states = [+3]

    class C:
        atomic_number = 6
        symbol = "C"
        name = "Carbon"

        category = "nonmetal"
        electron_configuration = "[He] 2s2 2p2"

        group = 14
        period = 2

        molar_mass = 12.01
        electronegativity = 2.55
        oxidation_states = [-4, +2, +4]

    class N:
        atomic_number = 7
        symbol = "N"
        name = "Nitrogen"

        category = "nonmetal"
        electron_configuration = "[He] 2s2 2p3"

        group = 15
        period = 2

        molar_mass = 14.01
        electronegativity = 3.04
        oxidation_states = [-3, +3, +5]

    class O:
        atomic_number = 8
        symbol = "O"
        name = "Oxygen"

        category = "nonmetal"
        electron_configuration = "[He] 2s2 2p4"

        group = 16
        period = 2

        molar_mass = 16.00
        electronegativity = 3.44
        oxidation_states = [-2]

    class F:
        atomic_number = 9
        symbol = "F"
        name = "Fluorine"

        category = "halogen"
        electron_configuration = "[He] 2s2 2p5"

        group = 17
        period = 2

        molar_mass = 19.00
        electronegativity = 3.98
        oxidation_states = [-1]

    class Ne:
        atomic_number = 10
        symbol = "Ne"
        name = "Neon"

        category = "noble gas"
        electron_configuration = "[He] 2s2 2p6"

        group = 18
        period = 2

        molar_mass = 20.18
        electronegativity = None
        oxidation_states = [0]

    class Na:
        atomic_number = 11
        symbol = "Na"
        name = "Sodium"

        category = "alkali metal"
        electron_configuration = "[Ne] 3s1"

        group = 1
        period = 3

        molar_mass = 22.99
        electronegativity = 0.93
        oxidation_states = [+1]

    class Mg:
        atomic_number = 12
        symbol = "Mg"
        name = "Magnesium"

        category = "alkaline earth metal"
        electron_configuration = "[Ne] 3s2"

        group = 2
        period = 3

        molar_mass = 24.31
        electronegativity = 1.31
        oxidation_states = [+2]

    class Al:
        atomic_number = 13
        symbol = "Al"
        name = "Aluminium"

        category = "post-transition metal"
        electron_configuration = "[Ne] 3s2 3p1"

        group = 13
        period = 3

        molar_mass = 26.98
        electronegativity = 1.61
        oxidation_states = [+3]

    class Si:
        atomic_number = 14
        symbol = "Si"
        name = "Silicon"

        category = "metalloid"
        electron_configuration = "[Ne] 3s2 3p2"

        group = 14
        period = 3

        molar_mass = 28.09
        electronegativity = 1.90
        oxidation_states = [-4, +2, +4]

    class P:
        atomic_number = 15
        symbol = "P"
        name = "Phosphorus"

        category = "nonmetal"
        electron_configuration = "[Ne] 3s2 3p3"

        group = 15
        period = 3

        molar_mass = 30.97
        electronegativity = 2.19
        oxidation_states = [-3, +3, +5]

    class S:
        atomic_number = 16
        symbol = "S"
        name = "Sulfur"

        category = "nonmetal"
        electron_configuration = "[Ne] 3s2 3p4"

        group = 16
        period = 3

        molar_mass = 32.06
        electronegativity = 2.58
        oxidation_states = [-2, +4, +6]

    class Cl:
        atomic_number = 17
        symbol = "Cl"
        name = "Chlorine"

        category = "halogen"
        electron_configuration = "[Ne] 3s2 3p5"

        group = 17
        period = 3

        molar_mass = 35.45
        electronegativity = 3.16
        oxidation_states = [-1, +1, +3, +5, +7]

    class Ar:
        atomic_number = 18
        symbol = "Ar"
        name = "Argon"

        category = "noble gas"
        electron_configuration = "[Ne] 3s2 3p6"

        group = 18
        period = 3

        molar_mass = 39.95
        electronegativity = None
        oxidation_states = [0]

    class K:
        atomic_number = 19
        symbol = "K"
        name = "Potassium"

        category = "alkali metal"
        electron_configuration = "[Ar] 4s1"

        group = 1
        period = 4

        molar_mass = 39.10
        electronegativity = 0.82
        oxidation_states = [+1]

    class Ca:
        atomic_number = 20
        symbol = "Ca"
        name = "Calcium"

        category = "alkaline earth metal"
        electron_configuration = "[Ar] 4s2"

        group = 2
        period = 4

        molar_mass = 40.08
        electronegativity = 1.00
        oxidation_states = [+2]

    class Sc:
        atomic_number = 21
        symbol = "Sc"
        name = "Scandium"

        category = "transition metal"
        electron_configuration = "[Ar] 3d1 4s2"

        group = 3
        period = 4

        molar_mass = 44.96
        electronegativity = 1.36
        oxidation_states = [+3]

    class Ti:
        atomic_number = 22
        symbol = "Ti"
        name = "Titanium"

        category = "transition metal"
        electron_configuration = "[Ar] 3d2 4s2"

        group = 4
        period = 4

        molar_mass = 47.87
        electronegativity = 1.54
        oxidation_states = [+2, +3, +4]

    class V:
        atomic_number = 23
        symbol = "V"
        name = "Vanadium"

        category = "transition metal"
        electron_configuration = "[Ar] 3d3 4s2"

        group = 5
        period = 4

        molar_mass = 50.94
        electronegativity = 1.63
        oxidation_states = [+2, +3, +4, +5]

    class Cr:
        atomic_number = 24
        symbol = "Cr"
        name = "Chromium"

        category = "transition metal"
        electron_configuration = "[Ar] 3d5 4s1"

        group = 6
        period = 4

        molar_mass = 52.00
        electronegativity = 1.66
        oxidation_states = [+2, +3, +6]

    class Mn:
        atomic_number = 25
        symbol = "Mn"
        name = "Manganese"

        category = "transition metal"
        electron_configuration = "[Ar] 3d5 4s2"

        group = 7
        period = 4

        molar_mass = 54.94
        electronegativity = 1.55
        oxidation_states = [+2, +3, +4, +6, +7]

    class Fe:
        atomic_number = 26
        symbol = "Fe"
        name = "Iron"

        category = "transition metal"
        electron_configuration = "[Ar] 3d6 4s2"

        group = 8
        period = 4

        molar_mass = 55.85
        electronegativity = 1.83
        oxidation_states = [+2, +3]

    class Co:
        atomic_number = 27
        symbol = "Co"
        name = "Cobalt"

        category = "transition metal"
        electron_configuration = "[Ar] 3d7 4s2"

        group = 9
        period = 4

        molar_mass = 58.93
        electronegativity = 1.88
        oxidation_states = [+2, +3]

    class Ni:
        atomic_number = 28
        symbol = "Ni"
        name = "Nickel"

        category = "transition metal"
        electron_configuration = "[Ar] 3d8 4s2"

        group = 10
        period = 4

        molar_mass = 58.69
        electronegativity = 1.91
        oxidation_states = [+2, +3]

    class Cu:
        atomic_number = 29
        symbol = "Cu"
        name = "Copper"

        category = "transition metal"
        electron_configuration = "[Ar] 3d10 4s1"

        group = 11
        period = 4

        molar_mass = 63.55
        electronegativity = 1.90
        oxidation_states = [+1, +2]

    class Zn:
        atomic_number = 30
        symbol = "Zn"
        name = "Zinc"

        category = "transition metal"
        electron_configuration = "[Ar] 3d10 4s2"

        group = 12
        period = 4

        molar_mass = 65.38
        electronegativity = 1.65
        oxidation_states = [+2]

    class Ga:
        atomic_number = 31
        symbol = "Ga"
        name = "Gallium"

        category = "post-transition metal"
        electron_configuration = "[Ar] 3d10 4s2 4p1"

        group = 13
        period = 4

        molar_mass = 69.72
        electronegativity = 1.81
        oxidation_states = [+3]

    class Ge:
        atomic_number = 32
        symbol = "Ge"
        name = "Germanium"

        category = "metalloid"
        electron_configuration = "[Ar] 3d10 4s2 4p2"

        group = 14
        period = 4

        molar_mass = 72.63
        electronegativity = 2.01
        oxidation_states = [+2, +4]

    class As:
        atomic_number = 33
        symbol = "As"
        name = "Arsenic"

        category = "metalloid"
        electron_configuration = "[Ar] 3d10 4s2 4p3"

        group = 15
        period = 4

        molar_mass = 74.92
        electronegativity = 2.18
        oxidation_states = [-3, +3, +5]

    class Se:
        atomic_number = 34
        symbol = "Se"
        name = "Selenium"

        category = "nonmetal"
        electron_configuration = "[Ar] 3d10 4s2 4p4"

        group = 16
        period = 4

        molar_mass = 78.97
        electronegativity = 2.55
        oxidation_states = [-2, +4, +6]

    class Br:
        atomic_number = 35
        symbol = "Br"
        name = "Bromine"

        category = "halogen"
        electron_configuration = "[Ar] 3d10 4s2 4p5"

        group = 17
        period = 4

        molar_mass = 79.90
        electronegativity = 2.96
        oxidation_states = [-1, +1, +3, +5]

    class Kr:
        atomic_number = 36
        symbol = "Kr"
        name = "Krypton"

        category = "noble gas"
        electron_configuration = "[Ar] 3d10 4s2 4p6"

        group = 18
        period = 4

        molar_mass = 83.80
        electronegativity = 3.00
        oxidation_states = [0]

    class Rb:
        atomic_number = 37
        symbol = "Rb"
        name = "Rubidium"

        category = "alkali metal"
        electron_configuration = "[Kr] 5s1"

        group = 1
        period = 5

        molar_mass = 85.47
        electronegativity = 0.82
        oxidation_states = [+1]

    class Sr:
        atomic_number = 38
        symbol = "Sr"
        name = "Strontium"

        category = "alkaline earth metal"
        electron_configuration = "[Kr] 5s2"

        group = 2
        period = 5

        molar_mass = 87.62
        electronegativity = 0.95
        oxidation_states = [+2]

    class Y:
        atomic_number = 39
        symbol = "Y"
        name = "Yttrium"

        category = "transition metal"
        electron_configuration = "[Kr] 4d1 5s2"

        group = 3
        period = 5

        molar_mass = 88.91
        electronegativity = 1.22
        oxidation_states = [+3]

    class Zr:
        atomic_number = 40
        symbol = "Zr"
        name = "Zirconium"

        category = "transition metal"
        electron_configuration = "[Kr] 4d2 5s2"

        group = 4
        period = 5

        molar_mass = 91.22
        electronegativity = 1.33
        oxidation_states = [+4]

    class Nb:
        atomic_number = 41
        symbol = "Nb"
        name = "Niobium"

        category = "transition metal"
        electron_configuration = "[Kr] 4d4 5s1"

        group = 5
        period = 5

        molar_mass = 92.91
        electronegativity = 1.6
        oxidation_states = [+3, +5]

    class Mo:
        atomic_number = 42
        symbol = "Mo"
        name = "Molybdenum"

        category = "transition metal"
        electron_configuration = "[Kr] 4d5 5s1"

        group = 6
        period = 5

        molar_mass = 95.95
        electronegativity = 2.16
        oxidation_states = [+4, +6]

    class Tc:
        atomic_number = 43
        symbol = "Tc"
        name = "Technetium"

        category = "transition metal"
        electron_configuration = "[Kr] 4d5 5s2"

        group = 7
        period = 5

        molar_mass = 98.00
        electronegativity = 1.9
        oxidation_states = [+4, +7]

    class Ru:
        atomic_number = 44
        symbol = "Ru"
        name = "Ruthenium"

        category = "transition metal"
        electron_configuration = "[Kr] 4d7 5s1"

        group = 8
        period = 5

        molar_mass = 101.07
        electronegativity = 2.2
        oxidation_states = [+3, +4, +8]

    class Rh:
        atomic_number = 45
        symbol = "Rh"
        name = "Rhodium"

        category = "transition metal"
        electron_configuration = "[Kr] 4d8 5s1"

        group = 9
        period = 5

        molar_mass = 102.91
        electronegativity = 2.28
        oxidation_states = [+3]

    class Pd:
        atomic_number = 46
        symbol = "Pd"
        name = "Palladium"

        category = "transition metal"
        electron_configuration = "[Kr] 4d10"

        group = 10
        period = 5

        molar_mass = 106.42
        electronegativity = 2.20
        oxidation_states = [+2, +4]

    class Ag:
        atomic_number = 47
        symbol = "Ag"
        name = "Silver"

        category = "transition metal"
        electron_configuration = "[Kr] 4d10 5s1"

        group = 11
        period = 5

        molar_mass = 107.87
        electronegativity = 1.93
        oxidation_states = [+1]

    class Cd:
        atomic_number = 48
        symbol = "Cd"
        name = "Cadmium"

        category = "transition metal"
        electron_configuration = "[Kr] 4d10 5s2"

        group = 12
        period = 5

        molar_mass = 112.41
        electronegativity = 1.69
        oxidation_states = [+2]

    class In:
        atomic_number = 49
        symbol = "In"
        name = "Indium"

        category = "post-transition metal"
        electron_configuration = "[Kr] 4d10 5s2 5p1"

        group = 13
        period = 5

        molar_mass = 114.82
        electronegativity = 1.78
        oxidation_states = [+3]

    class Sn:
        atomic_number = 50
        symbol = "Sn"
        name = "Tin"

        category = "post-transition metal"
        electron_configuration = "[Kr] 4d10 5s2 5p2"

        group = 14
        period = 5

        molar_mass = 118.71
        electronegativity = 1.96
        oxidation_states = [+2, +4]

    class Sb:
        atomic_number = 51
        symbol = "Sb"
        name = "Antimony"

        category = "metalloid"
        electron_configuration = "[Kr] 4d10 5s2 5p3"

        group = 15
        period = 5

        molar_mass = 121.76
        electronegativity = 2.05
        oxidation_states = [-3, +3, +5]

    class Te:
        atomic_number = 52
        symbol = "Te"
        name = "Tellurium"

        category = "metalloid"
        electron_configuration = "[Kr] 4d10 5s2 5p4"

        group = 16
        period = 5

        molar_mass = 127.60
        electronegativity = 2.1
        oxidation_states = [-2, +4, +6]

    class I:
        atomic_number = 53
        symbol = "I"
        name = "Iodine"

        category = "halogen"
        electron_configuration = "[Kr] 4d10 5s2 5p5"

        group = 17
        period = 5

        molar_mass = 126.90
        electronegativity = 2.66
        oxidation_states = [-1, +1, +3, +5, +7]

    class Xe:
        atomic_number = 54
        symbol = "Xe"
        name = "Xenon"

        category = "noble gas"
        electron_configuration = "[Kr] 4d10 5s2 5p6"

        group = 18
        period = 5

        molar_mass = 131.29
        electronegativity = 2.6
        oxidation_states = [0, +2, +4, +6, +8]

    class Cs:
        atomic_number = 55
        symbol = "Cs"
        name = "Caesium"

        category = "alkali metal"
        electron_configuration = "[Xe] 6s1"

        group = 1
        period = 6

        molar_mass = 132.91
        electronegativity = 0.79
        oxidation_states = [+1]

    class Ba:
        atomic_number = 56
        symbol = "Ba"
        name = "Barium"

        category = "alkaline earth metal"
        electron_configuration = "[Xe] 6s2"

        group = 2
        period = 6

        molar_mass = 137.33
        electronegativity = 0.89
        oxidation_states = [+2]

    class La:
        atomic_number = 57
        symbol = "La"
        name = "Lanthanum"

        category = "lanthanide"
        electron_configuration = "[Xe] 5d1 6s2"

        group = 3
        period = 6

        molar_mass = 138.91
        electronegativity = 1.10
        oxidation_states = [+3]

    class Ce:
        atomic_number = 58
        symbol = "Ce"
        name = "Cerium"

        category = "lanthanide"
        electron_configuration = "[Xe] 4f1 5d1 6s2"

        group = None
        period = 6

        molar_mass = 140.12
        electronegativity = 1.12
        oxidation_states = [+3, +4]

    class Pr:
        atomic_number = 59
        symbol = "Pr"
        name = "Praseodymium"

        category = "lanthanide"
        electron_configuration = "[Xe] 4f3 6s2"

        group = None
        period = 6

        molar_mass = 140.91
        electronegativity = 1.13
        oxidation_states = [+3, +4]

    class Nd:
        atomic_number = 60
        symbol = "Nd"
        name = "Neodymium"

        category = "lanthanide"
        electron_configuration = "[Xe] 4f4 6s2"

        group = None
        period = 6

        molar_mass = 144.24
        electronegativity = 1.14
        oxidation_states = [+3]

    class Pm:
        atomic_number = 61
        symbol = "Pm"
        name = "Promethium"

        category = "lanthanide"
        electron_configuration = "[Xe] 4f5 6s2"

        group = None
        period = 6

        molar_mass = 145.00
        electronegativity = 1.13
        oxidation_states = [+3]

    class Sm:
        atomic_number = 62
        symbol = "Sm"
        name = "Samarium"

        category = "lanthanide"
        electron_configuration = "[Xe] 4f6 6s2"

        group = None
        period = 6

        molar_mass = 150.36
        electronegativity = 1.17
        oxidation_states = [+2, +3]

    class Eu:
        atomic_number = 63
        symbol = "Eu"
        name = "Europium"

        category = "lanthanide"
        electron_configuration = "[Xe] 4f7 6s2"

        group = None
        period = 6

        molar_mass = 151.96
        electronegativity = 1.20
        oxidation_states = [+2, +3]

    class Gd:
        atomic_number = 64
        symbol = "Gd"
        name = "Gadolinium"

        category = "lanthanide"
        electron_configuration = "[Xe] 4f7 5d1 6s2"

        group = None
        period = 6

        molar_mass = 157.25
        electronegativity = 1.20
        oxidation_states = [+3]

    class Tb:
        atomic_number = 65
        symbol = "Tb"
        name = "Terbium"

        category = "lanthanide"
        electron_configuration = "[Xe] 4f9 6s2"

        group = None
        period = 6

        molar_mass = 158.93
        electronegativity = 1.10
        oxidation_states = [+3, +4]

    class Dy:
        atomic_number = 66
        symbol = "Dy"
        name = "Dysprosium"

        category = "lanthanide"
        electron_configuration = "[Xe] 4f10 6s2"

        group = None
        period = 6

        molar_mass = 162.50
        electronegativity = 1.22
        oxidation_states = [+3]

    class Ho:
        atomic_number = 67
        symbol = "Ho"
        name = "Holmium"

        category = "lanthanide"
        electron_configuration = "[Xe] 4f11 6s2"

        group = None
        period = 6

        molar_mass = 164.93
        electronegativity = 1.23
        oxidation_states = [+3]

    class Er:
        atomic_number = 68
        symbol = "Er"
        name = "Erbium"

        category = "lanthanide"
        electron_configuration = "[Xe] 4f12 6s2"

        group = None
        period = 6

        molar_mass = 167.26
        electronegativity = 1.24
        oxidation_states = [+3]

    class Tm:
        atomic_number = 69
        symbol = "Tm"
        name = "Thulium"

        category = "lanthanide"
        electron_configuration = "[Xe] 4f13 6s2"

        group = None
        period = 6

        molar_mass = 168.93
        electronegativity = 1.25
        oxidation_states = [+3]

    class Yb:
        atomic_number = 70
        symbol = "Yb"
        name = "Ytterbium"

        category = "lanthanide"
        electron_configuration = "[Xe] 4f14 6s2"

        group = None
        period = 6

        molar_mass = 173.05
        electronegativity = 1.10
        oxidation_states = [+2, +3]

    class Lu:
        atomic_number = 71
        symbol = "Lu"
        name = "Lutetium"

        category = "lanthanide"
        electron_configuration = "[Xe] 4f14 5d1 6s2"

        group = 3
        period = 6

        molar_mass = 174.97
        electronegativity = 1.27
        oxidation_states = [+3]

    class Hf:
        atomic_number = 72
        symbol = "Hf"
        name = "Hafnium"

        category = "transition metal"
        electron_configuration = "[Xe] 4f14 5d2 6s2"

        group = 4
        period = 6

        molar_mass = 178.49
        electronegativity = 1.3
        oxidation_states = [+4]

    class Ta:
        atomic_number = 73
        symbol = "Ta"
        name = "Tantalum"

        category = "transition metal"
        electron_configuration = "[Xe] 4f14 5d3 6s2"

        group = 5
        period = 6

        molar_mass = 180.95
        electronegativity = 1.5
        oxidation_states = [+5]

    class W:
        atomic_number = 74
        symbol = "W"
        name = "Tungsten"

        category = "transition metal"
        electron_configuration = "[Xe] 4f14 5d4 6s2"

        group = 6
        period = 6

        molar_mass = 183.84
        electronegativity = 2.36
        oxidation_states = [+6]

    class Re:
        atomic_number = 75
        symbol = "Re"
        name = "Rhenium"

        category = "transition metal"
        electron_configuration = "[Xe] 4f14 5d5 6s2"

        group = 7
        period = 6

        molar_mass = 186.21
        electronegativity = 1.9
        oxidation_states = [+4, +7]

    class Os:
        atomic_number = 76
        symbol = "Os"
        name = "Osmium"

        category = "transition metal"
        electron_configuration = "[Xe] 4f14 5d6 6s2"

        group = 8
        period = 6

        molar_mass = 190.23
        electronegativity = 2.2
        oxidation_states = [+4, +6, +8]

    class Ir:
        atomic_number = 77
        symbol = "Ir"
        name = "Iridium"

        category = "transition metal"
        electron_configuration = "[Xe] 4f14 5d7 6s2"

        group = 9
        period = 6

        molar_mass = 192.22
        electronegativity = 2.20
        oxidation_states = [+3, +4, +6]

    class Pt:
        atomic_number = 78
        symbol = "Pt"
        name = "Platinum"

        category = "transition metal"
        electron_configuration = "[Xe] 4f14 5d9 6s1"

        group = 10
        period = 6

        molar_mass = 195.08
        electronegativity = 2.28
        oxidation_states = [+2, +4]

    class Au:
        atomic_number = 79
        symbol = "Au"
        name = "Gold"

        category = "transition metal"
        electron_configuration = "[Xe] 4f14 5d10 6s1"

        group = 11
        period = 6

        molar_mass = 196.97
        electronegativity = 2.54
        oxidation_states = [+1, +3]

    class Hg:
        atomic_number = 80
        symbol = "Hg"
        name = "Mercury"

        category = "transition metal"
        electron_configuration = "[Xe] 4f14 5d10 6s2"

        group = 12
        period = 6

        molar_mass = 200.59
        electronegativity = 2.00
        oxidation_states = [+1, +2]

    class Tl:
        atomic_number = 81
        symbol = "Tl"
        name = "Thallium"

        category = "post-transition metal"
        electron_configuration = "[Xe] 4f14 5d10 6s2 6p1"

        group = 13
        period = 6

        molar_mass = 204.38
        electronegativity = 1.62
        oxidation_states = [+1, +3]

    class Pb:
        atomic_number = 82
        symbol = "Pb"
        name = "Lead"

        category = "post-transition metal"
        electron_configuration = "[Xe] 4f14 5d10 6s2 6p2"

        group = 14
        period = 6

        molar_mass = 207.2
        electronegativity = 2.33
        oxidation_states = [+2, +4]

    class Bi:
        atomic_number = 83
        symbol = "Bi"
        name = "Bismuth"

        category = "post-transition metal"
        electron_configuration = "[Xe] 4f14 5d10 6s2 6p3"

        group = 15
        period = 6

        molar_mass = 208.98
        electronegativity = 2.02
        oxidation_states = [+3, +5]

    class Po:
        atomic_number = 84
        symbol = "Po"
        name = "Polonium"

        category = "metalloid"
        electron_configuration = "[Xe] 4f14 5d10 6s2 6p4"

        group = 16
        period = 6

        molar_mass = 209.00
        electronegativity = 2.0
        oxidation_states = [-2, +2, +4, +6]

    class At:
        atomic_number = 85
        symbol = "At"
        name = "Astatine"

        category = "halogen"
        electron_configuration = "[Xe] 4f14 5d10 6s2 6p5"

        group = 17
        period = 6

        molar_mass = 210.00
        electronegativity = 2.2
        oxidation_states = [-1, +1, +3, +5, +7]

    class Rn:
        atomic_number = 86
        symbol = "Rn"
        name = "Radon"

        category = "noble gas"
        electron_configuration = "[Xe] 4f14 5d10 6s2 6p6"

        group = 18
        period = 6

        molar_mass = 222.00
        electronegativity = None
        oxidation_states = [0]

    class Fr:
        atomic_number = 87
        symbol = "Fr"
        name = "Francium"

        category = "alkali metal"
        electron_configuration = "[Rn] 7s1"

        group = 1
        period = 7

        molar_mass = 223
        electronegativity = 0.7
        oxidation_states = [+1]

    class Ra:
        atomic_number = 88
        symbol = "Ra"
        name = "Radium"

        category = "alkaline earth metal"
        electron_configuration = "[Rn] 7s2"

        group = 2
        period = 7

        molar_mass = 226
        electronegativity = 0.9
        oxidation_states = [+2]

    class Ac:
        atomic_number = 89
        symbol = "Ac"
        name = "Actinium"

        category = "actinide"
        electron_configuration = "[Rn] 6d1 7s2"

        group = 3
        period = 7

        molar_mass = 227
        electronegativity = 1.1
        oxidation_states = [+3]

    class Th:
        atomic_number = 90
        symbol = "Th"
        name = "Thorium"

        category = "actinide"
        electron_configuration = "[Rn] 6d2 7s2"

        group = None
        period = 7

        molar_mass = 232.04
        electronegativity = 1.3
        oxidation_states = [+4]

    class Pa:
        atomic_number = 91
        symbol = "Pa"
        name = "Protactinium"

        category = "actinide"
        electron_configuration = "[Rn] 5f2 6d1 7s2"

        group = None
        period = 7

        molar_mass = 231.04
        electronegativity = 1.5
        oxidation_states = [+5]

    class U:
        atomic_number = 92
        symbol = "U"
        name = "Uranium"

        category = "actinide"
        electron_configuration = "[Rn] 5f3 6d1 7s2"

        group = None
        period = 7

        molar_mass = 238.03
        electronegativity = 1.38
        oxidation_states = [+3, +4, +5, +6]

    class Np:
        atomic_number = 93
        symbol = "Np"
        name = "Neptunium"

        category = "actinide"
        electron_configuration = "[Rn] 5f4 6d1 7s2"

        group = None
        period = 7

        molar_mass = 237
        electronegativity = 1.36
        oxidation_states = [+3, +4, +5, +6, +7]

    class Pu:
        atomic_number = 94
        symbol = "Pu"
        name = "Plutonium"

        category = "actinide"
        electron_configuration = "[Rn] 5f6 7s2"

        group = None
        period = 7

        molar_mass = 244
        electronegativity = 1.28
        oxidation_states = [+3, +4, +5, +6]

    class Am:
        atomic_number = 95
        symbol = "Am"
        name = "Americium"

        category = "actinide"
        electron_configuration = "[Rn] 5f7 7s2"

        group = None
        period = 7

        molar_mass = 243
        electronegativity = 1.3
        oxidation_states = [+3, +4, +6]

    class Cm:
        atomic_number = 96
        symbol = "Cm"
        name = "Curium"

        category = "actinide"
        electron_configuration = "[Rn] 5f7 6d1 7s2"

        group = None
        period = 7

        molar_mass = 247
        electronegativity = 1.3
        oxidation_states = [+3, +4]

    class Bk:
        atomic_number = 97
        symbol = "Bk"
        name = "Berkelium"

        category = "actinide"
        electron_configuration = "[Rn] 5f9 7s2"

        group = None
        period = 7

        molar_mass = 247
        electronegativity = 1.3
        oxidation_states = [+3, +4]

    class Cf:
        atomic_number = 98
        symbol = "Cf"
        name = "Californium"

        category = "actinide"
        electron_configuration = "[Rn] 5f10 7s2"

        group = None
        period = 7

        molar_mass = 251
        electronegativity = 1.3
        oxidation_states = [+3, +4]

    class Es:
        atomic_number = 99
        symbol = "Es"
        name = "Einsteinium"

        category = "actinide"
        electron_configuration = "[Rn] 5f11 7s2"

        group = None
        period = 7

        molar_mass = 252
        electronegativity = 1.3
        oxidation_states = [+3]

    class Fm:
        atomic_number = 100
        symbol = "Fm"
        name = "Fermium"

        category = "actinide"
        electron_configuration = "[Rn] 5f12 7s2"

        group = None
        period = 7

        molar_mass = 257
        electronegativity = 1.3
        oxidation_states = [+3]

    class Md:
        atomic_number = 101
        symbol = "Md"
        name = "Mendelevium"

        category = "actinide"
        electron_configuration = "[Rn] 5f13 7s2"

        group = None
        period = 7

        molar_mass = 258
        electronegativity = 1.3
        oxidation_states = [+2, +3]

    class No:
        atomic_number = 102
        symbol = "No"
        name = "Nobelium"

        category = "actinide"
        electron_configuration = "[Rn] 5f14 7s2"

        group = None
        period = 7

        molar_mass = 259
        electronegativity = 1.3
        oxidation_states = [+2, +3]

    class Lr:
        atomic_number = 103
        symbol = "Lr"
        name = "Lawrencium"

        category = "actinide"
        electron_configuration = "[Rn] 5f14 7s2 7p1"

        group = 3
        period = 7

        molar_mass = 262
        electronegativity = 1.3
        oxidation_states = [+3]

    class Rf:
        atomic_number = 104
        symbol = "Rf"
        name = "Rutherfordium"

        category = "transition metal"
        electron_configuration = "[Rn] 5f14 6d2 7s2"

        group = 4
        period = 7

        molar_mass = 267
        electronegativity = None
        oxidation_states = [+4]

    class Db:
        atomic_number = 105
        symbol = "Db"
        name = "Dubnium"

        category = "transition metal"
        electron_configuration = "[Rn] 5f14 6d3 7s2"

        group = 5
        period = 7

        molar_mass = 270
        electronegativity = None
        oxidation_states = [+5]

    class Sg:
        atomic_number = 106
        symbol = "Sg"
        name = "Seaborgium"

        category = "transition metal"
        electron_configuration = "[Rn] 5f14 6d4 7s2"

        group = 6
        period = 7

        molar_mass = 271
        electronegativity = None
        oxidation_states = [+6]

    class Bh:
        atomic_number = 107
        symbol = "Bh"
        name = "Bohrium"

        category = "transition metal"
        electron_configuration = "[Rn] 5f14 6d5 7s2"

        group = 7
        period = 7

        molar_mass = 270
        electronegativity = None
        oxidation_states = [+7]

    class Hs:
        atomic_number = 108
        symbol = "Hs"
        name = "Hassium"

        category = "transition metal"
        electron_configuration = "[Rn] 5f14 6d6 7s2"

        group = 8
        period = 7

        molar_mass = 277
        electronegativity = None
        oxidation_states = [+8]

    class Mt:
        atomic_number = 109
        symbol = "Mt"
        name = "Meitnerium"

        category = "transition metal"
        electron_configuration = "[Rn] 5f14 6d7 7s2"

        group = 9
        period = 7

        molar_mass = 278
        electronegativity = None
        oxidation_states = [+1, +3]

    class Ds:
        atomic_number = 110
        symbol = "Ds"
        name = "Darmstadtium"

        category = "transition metal"
        electron_configuration = "[Rn] 5f14 6d8 7s2"

        group = 10
        period = 7

        molar_mass = 281
        electronegativity = None
        oxidation_states = [+2, +4]

    class Rg:
        atomic_number = 111
        symbol = "Rg"
        name = "Roentgenium"

        category = "transition metal"
        electron_configuration = "[Rn] 5f14 6d9 7s2"

        group = 11
        period = 7

        molar_mass = 282
        electronegativity = None
        oxidation_states = [+1, +3]

    class Cn:
        atomic_number = 112
        symbol = "Cn"
        name = "Copernicium"

        category = "transition metal"
        electron_configuration = "[Rn] 5f14 6d10 7s2"

        group = 12
        period = 7

        molar_mass = 285
        electronegativity = None
        oxidation_states = [+2, +4]

    class Nh:
        atomic_number = 113
        symbol = "Nh"
        name = "Nihonium"

        category = "post-transition metal"
        electron_configuration = "[Rn] 5f14 6d10 7s2 7p1"

        group = 13
        period = 7

        molar_mass = 286
        electronegativity = None
        oxidation_states = [+1, +3]

    class Fl:
        atomic_number = 114
        symbol = "Fl"
        name = "Flerovium"

        category = "post-transition metal"
        electron_configuration = "[Rn] 5f14 6d10 7s2 7p2"

        group = 14
        period = 7

        molar_mass = 289
        electronegativity = None
        oxidation_states = [+2, +4]

    class Mc:
        atomic_number = 115
        symbol = "Mc"
        name = "Moscovium"

        category = "post-transition metal"
        electron_configuration = "[Rn] 5f14 6d10 7s2 7p3"

        group = 15
        period = 7

        molar_mass = 288
        electronegativity = None
        oxidation_states = [+1, +3, +5]

    class Lv:
        atomic_number = 116
        symbol = "Lv"
        name = "Livermorium"

        category = "post-transition metal"
        electron_configuration = "[Rn] 5f14 6d10 7s2 7p4"

        group = 16
        period = 7

        molar_mass = 293
        electronegativity = None
        oxidation_states = [+2, +4, +6]

    class Ts:
        atomic_number = 117
        symbol = "Ts"
        name = "Tennessine"

        category = "halogen"
        electron_configuration = "[Rn] 5f14 6d10 7s2 7p5"

        group = 17
        period = 7

        molar_mass = 294
        electronegativity = None
        oxidation_states = [-1, +1, +3, +5, +7]

    class Og:
        atomic_number = 118
        symbol = "Og"
        name = "Oganesson"

        category = "noble gas"
        electron_configuration = "[Rn] 5f14 6d10 7s2 7p6"

        group = 18
        period = 7

        molar_mass = 294
        electronegativity = None
        oxidation_states = [0]





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
