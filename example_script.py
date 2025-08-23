import chem
from chem import PeriodicTable as pt, Molekyle
import sympy


def using_the_periodic_table():
    He = pt.He()
    print(He.atomic_number)
    print(He.name)
    print(He.symbol)
    print(He.molar_mass)

    # ----
    print()
    # ----

    H2 = pt.H(2)
    print(H2.molar_mass)
    print(H2.total_molar_mass())

    # ----
    print()
    # ----

    NH3 = Molekyle(
        pt.N(),  # pt.N(1)
        pt.H(3)
    )
    NH3.state = 'g'
    print(NH3.get_thermal_properties())


def molar_mass():
    C2H5OH = chem.Molekyle(
        pt.C(2),
        pt.H(6),
        pt.O(1)
    )

    # Let's say that we have 2.3 grams of ethanol
    grams = 2.3
    C2H5OH.grams = grams
    print(f'{C2H5OH.grams = }')

    # We can then calculate the amount of mole ethanol
    # we have based on the molar mass of the molekyle.
    # The molar mass is already calculated, and can be accessed by writing:
    C2H5OH.total_molar_mass()
    print(f'{C2H5OH.total_molar_mass() = }')

    # And we can just get the mol by writing:
    C2H5OH.calculate_mol()
    print(f'{C2H5OH.mol = }')
    # This calculation can only be done if we have given the amount of grams

    # In the same way we can write:
    C2H5OH.grams = 0
    # So we calculate grams based on mol
    C2H5OH.calculate_gram()
    print(f'{C2H5OH.grams = }')
    # Which calculates the amount of grams based on the mol given
    # This should be the same as the amount of grams specified in the start

    # If you don't want the methods to be black boxes then I suggest
    # taking a look into the Molekyle class and Element class


def using_a_formula():
    V = 1    # L
    p = 5.4  # atm
    T = 315  # K

    IGL = chem.IdealGasLaw()
    IGL.V = V
    IGL.p = p
    IGL.T = T
    print(f'{IGL.compute() = }')

    IGL = chem.IdealGasLaw(V=V, T=T, p=p)
    print(f'{IGL.compute() = }')

    result = chem.IdealGasLaw(V=V, T=T, p=p).compute()
    print(f'{result = }')

    isolated = IGL.solve(IGL.n)
    print(f'{isolated = }')

    # ----
    print()
    # ----

    A = 0.2
    K = 4.6e-4

    Eq = chem.Equilibrium()
    x = sympy.symbols('x')

    Eq.K = K

    Eq.A = A - x
    Eq.B = 1

    Eq.C = x
    Eq.D = x

    result = Eq.compute()
    concentration = sympy.solve(result, x)

    print(f'{concentration = }')

    # ----
    print()
    # ----

    C = 25  # °C
    K = chem.Kelvin(C=C).compute()
    print(f'{K = }')

    K = 315
    C = chem.Kelvin(K=K).compute()
    print(f'{C = }')

    # ----
    print()
    # ----

    unit = chem.UnitConverter()
    kJ = 4.3  # kJ
    J = unit.convert(kJ, 'k', '')
    print(f'{J = }')

    volume = chem.VolumeConversion()
    L = 10
    m3 = volume.L_to_m3(L)
    print(f'{m3 = }')


if __name__ == '__main__':
    using_the_periodic_table()
