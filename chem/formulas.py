from .constants import Constants
from .formula_backbone import *
import sympy as sp


TYPE = Variable | int | float | sp.Symbol | None


class IdealGasLaw(Formula):
    """
    **Ideal Gas Law**::

        pV = nRT

        p : Pa          (pressure)
        V : m³          (volume)
        n : mol         (amount of substance)
        R : J/(mol·K)   (gas constant)
        T : K           (temperature)
    """
    p: TYPE
    V: TYPE
    n: TYPE
    R: TYPE
    T: TYPE

    def __init__(self, *, p=None, V=None, n=None, R=Constants.R_L, T=None):
        super().__init__(p, V, n, R, T)

        p, V, n, R, T = self.symbols

        self.formula = p*V - n*R*T


class CoulombsLaw(Formula):
    """
    **Coulombs Law**::

        F_e = k_e · |q1 · q2| / r²

        F_e : Force             (N)
        k_e : Coulomb constant  (N·m²/C²)
        q1  : Charge 1          (C)
        q2  : Charge 2          (C)
        r   : Distance          (m)
    """
    F_e: TYPE
    k_e: TYPE
    q1: TYPE
    q2: TYPE
    r: TYPE

    def __init__(self, *, F_e=None, k_e=Constants.k_e, q1=None, q2=None, r=None):
        super().__init__(F_e, k_e, q1, q2, r)

        F_e, k_e, q1, q2, r = self.symbols

        self.formula = F_e - k_e * sp.Abs(q1 * q2) / r**2


class Molarity(Formula):
    """
    **Molarity**::

        M = moles / V

        M     : Molarity    (mol/L)
        moles : Moles       (mol)
        V     : Volume      (L)
    """
    M: TYPE
    moles: TYPE
    V: TYPE

    def __init__(self, *, M=None, moles=None, V=None):
        super().__init__(M, moles, V)

        M, moles, V = self.symbols

        self.formula = M - moles / V


class Diluting(Formula):
    """
    **Diluting**::

        M_1·V_1 = M_2·V_2

        M : Molarity    (mol/L)
        V : Volume      (L)
    """
    M_1: TYPE
    V_1: TYPE
    M_2: TYPE
    V_2: TYPE

    def __init__(self, *, M_1=None, V_1=None, M_2=None, V_2=None):
        super().__init__(M_1, V_1, M_2, V_2)

        M_1, V_1, M_2, V_2 = self.symbols

        self.formula = M_1 * V_1 - M_2 * V_2


class BoylesLaw(Formula):
    """
    **Boyles Law**::

        P_1·V_1 = P_2·V_2

        P : Pressure (Pa)
        V : Volume   (M^3)
    """
    P_1: TYPE
    V_1: TYPE
    P_2: TYPE
    V_2: TYPE

    def __init__(self, *, P_1=None, V_1=None, P_2=None, V_2=None):
        super().__init__(P_1, V_1, P_2, V_2)

        P_1, V_1, P_2, V_2 = self.symbols

        self.formula = P_1 * V_1 - P_2 * V_2


class Kelvin(Formula):
    """
    **Kelvin**::

        K = C + 273.15

        K  : Kelvin
        C : Celsius
    """
    K: TYPE
    C: TYPE

    def __init__(self, *, K=None, C=None):
        super().__init__(K, C)

        K, C = self.symbols

        self.formula = K - (C + 273.15)


class CharlesLaw(Formula):
    """
    **Charles Law**::

        V_1 / T_1 = V_2 / T_2

        V : Volume      (M^3)
        T : Temperature (K)
    """
    V_1: TYPE
    T_1: TYPE
    V_2: TYPE
    T_2: TYPE

    def __init__(self, *, V_1=None, T_1=None, V_2=None, T_2=None):
        super().__init__(V_1, T_1, V_2, T_2)

        V_1, T_1, V_2, T_2 = self.symbols

        self.formula = V_1 / T_1 - V_2 / T_2


class CombinedGasLaw(Formula):
    """
    **Combined Gas Law**::

        P_1·V_1 / T_1 = P_2·V_2 / T_2

        P : Pressure (Pa)
        V : Volume (M^3)
        T : Temperature (K)
    """
    P_1: TYPE
    V_1: TYPE
    T_1: TYPE
    P_2: TYPE
    V_2: TYPE
    T_2: TYPE

    def __init__(self, *, P_1=None, V_1=None, T_1=None, P_2=None, V_2=None, T_2=None):
        super().__init__(P_1, V_1, T_1, P_2, V_2, T_2)

        P_1, V_1, T_1, P_2, V_2, T_2 = self.symbols

        self.formula = P_1 * V_1 / T_1 - P_2 * V_2 / T_2


class AvogadrosLaw(Formula):
    """
    **Avogadro's Law**::

        n_1 / V_1 = n_2 / V_2

        n : Mole   (mol)
        V : Volume (m^3)
    """
    n_1: TYPE
    V_1: TYPE
    n_2: TYPE
    V_2: TYPE

    def __init__(self, *, n_1=None, V_1=None, n_2=None, V_2=None):
        super().__init__(n_1, V_1, n_2, V_2)

        n_1, V_1, n_2, V_2 = self.symbols

        self.formula = n_1 / V_1 - n_2 / V_2


class KineticEnergy(Formula):
    """
    **Kinetic Energy**::

        E_k = 1/2 · m · v^2

        E_k : Kinetic Energy (J)
        m : Mass (kg)
        v : Speed (m/s)
    """
    E_k: TYPE
    m: TYPE
    v: TYPE

    def __init__(self, *, E_k=None, m=None, v=None):
        super().__init__(E_k, m, v)

        E_k, m, v = self.symbols

        self.formula = E_k - (1/2 * m * v**2)


class MostProbableSpeed(Formula):
    """
    **Most Probable Speed**::

        u_mp = sqrt(2·R·T/M)

        u_mp       : Most Probable Speed (m/s)
        R          : Ideal gas constant  (J/mol·K)
        T          : Temperature         (K)
        M : Molar mass          (kg/mol)
    """
    u_mp: TYPE
    R: TYPE
    T: TYPE
    molar_mass: TYPE

    def __init__(self, *, u_mp=None, R=Constants.R_J, T=None, M=None):
        super().__init__(u_mp, R, T, M)

        u_mp, R, T, M = self.symbols

        self.formula = u_mp - (sp.sqrt(2 * R * T / M))


class GrahamsLow(Formula):
    """
    **Grahams Low**::

        rate_2 / rate_1 = sqrt(M_1 / M_2)

        rate       : Effuse rate
        M : Molar mass (kg/mol)
    """
    rate_1: TYPE
    rate_2: TYPE
    molar_mass: TYPE

    def __init__(self, *, rate_1=None, rate_2=None, M_1=None, M_2=None):
        super().__init__(rate_1, rate_2, M_1, M_2)

        rate_1, rate_2, M_1, M_2 = self.symbols

        self.formula = rate_2 / rate_1 - sp.sqrt(M_1 / M_2)


class PressureConversion:
    atm_mmHg = 1/760
    atm_kPa = 1/101_325
    atm_bar = 1/1.013
    atm_psi = 1/14.70

    def mmHg_to_atm(self, mmHg):
        return mmHg * self.atm_mmHg

    def atm_to_mmHg(self, atm):
        return atm / self.atm_mmHg

    def kPa_to_atm(self, kPa):
        return kPa * self.atm_kPa

    def atm_to_kPa(self, atm):
        return atm / self.atm_kPa

    def bar_to_atm(self, bar):
        return bar * self.atm_bar

    def atm_to_bar(self, atm):
        return atm / self.atm_bar


class VolumeConversion:
    # Base unit: 1 L
    L_mL = 1000          # 1 L = 1000 mL
    L_cm3 = 1000         # 1 L = 1000 cm³
    L_dm3 = 1            # 1 L = 1 dm³
    L_m3 = 1/1000        # 1 L = 0.001 m³
    L_gal_US = 1/3.785   # 1 L ≈ 0.264 US gallons
    L_gal_UK = 1/4.546   # 1 L ≈ 0.220 UK gallons
    L_ft3 = 1/28.317     # 1 L ≈ 0.0353 ft³

    # --- to Liters ---
    def mL_to_L(self, mL): return mL / self.L_mL
    def cm3_to_L(self, cm3): return cm3 / self.L_cm3
    def dm3_to_L(self, dm3): return dm3 / self.L_dm3
    def m3_to_L(self, m3): return m3 * 1000
    def galUS_to_L(self, gal): return gal * 3.785
    def galUK_to_L(self, gal): return gal * 4.546
    def ft3_to_L(self, ft3): return ft3 * 28.317

    # --- from Liters ---
    def L_to_mL(self, L): return L * self.L_mL
    def L_to_cm3(self, L): return L * self.L_cm3
    def L_to_dm3(self, L): return L * self.L_dm3
    def L_to_m3(self, L): return L / 1000
    def L_to_galUS(self, L): return L * self.L_gal_US
    def L_to_galUK(self, L): return L * self.L_gal_UK
    def L_to_ft3(self, L): return L * self.L_ft3


class LengthConversion:
    # Base unit: meter (m)
    m_cm = 100            # 1 m = 100 cm
    m_mm = 1000           # 1 m = 1000 mm
    m_um = 1e6            # 1 m = 1,000,000 µm
    m_nm = 1e9            # 1 m = 1,000,000,000 nm
    m_pm = 1e12           # 1 m = 1,000,000,000,000 pm
    m_km = 1/1000         # 1 m = 0.001 km
    m_inch = 39.3701      # 1 m ≈ 39.37 in
    m_ft = 3.28084        # 1 m ≈ 3.28 ft
    m_yd = 1.09361        # 1 m ≈ 1.09 yd
    m_mile = 1/1609.34    # 1 m ≈ 0.000621 mi

    # --- to meters ---
    def cm_to_m(self, cm): return cm / self.m_cm
    def mm_to_m(self, mm): return mm / self.m_mm
    def um_to_m(self, um): return um / self.m_um
    def nm_to_m(self, nm): return nm / self.m_nm
    def pm_to_m(self, pm): return pm / self.m_pm
    def km_to_m(self, km): return km * 1000
    def inch_to_m(self, inch): return inch / self.m_inch
    def ft_to_m(self, ft): return ft / self.m_ft
    def yd_to_m(self, yd): return yd / self.m_yd
    def mile_to_m(self, mile): return mile * 1609.34

    # --- from meters ---
    def m_to_cm(self, m): return m * self.m_cm
    def m_to_mm(self, m): return m * self.m_mm
    def m_to_um(self, m): return m * self.m_um
    def m_to_nm(self, m): return m * self.m_nm
    def m_to_pm(self, m): return m * self.m_pm
    def m_to_km(self, m): return m * self.m_km
    def m_to_inch(self, m): return m * self.m_inch
    def m_to_ft(self, m): return m * self.m_ft
    def m_to_yd(self, m): return m * self.m_yd
    def m_to_mile(self, m): return m * self.m_mile


class UnitConverter:
    # SI prefix factors relative to base unit (1)
    PREFIXES = {
        "G": 1e9,      # giga
        "M": 1e6,      # mega
        "k": 1e3,      # kilo
        "h": 1e2,      # hecto
        "da": 1e1,     # deka
        "": 1,         # base unit
        "d": 1e-1,     # deci
        "c": 1e-2,     # centi
        "m": 1e-3,     # milli
        "u": 1e-6,     # micro (µ, but using 'u' for code safety)
        "n": 1e-9,     # nano
        "Å": 1e-10,    # Ångström
        "p": 1e-12     # pico
    }

    @classmethod
    def convert(cls, value: float, from_prefix: str, to_prefix: str) -> float:
        """
        Convert a value between two SI-prefixed units.

        Example: convert(5, "k", "m")  # 5 km -> 5,000,000 mm
        """
        if from_prefix not in cls.PREFIXES:
            raise ValueError(f"Invalid 'from' prefix: {from_prefix}")
        if to_prefix not in cls.PREFIXES:
            raise ValueError(f"Invalid 'to' prefix: {to_prefix}")

        # Convert to base unit first
        base_value = value * cls.PREFIXES[from_prefix]

        # Then convert to target prefix
        return base_value / cls.PREFIXES[to_prefix]


class EnergyConversion:
    # Constants
    J_per_eV = 1.602176634e-19   # 1 eV = 1.602e-19 J
    J_per_cal = 4.184            # 1 cal = 4.184 J
    J_per_kcal = 4184            # 1 kcal = 4184 J
    J_per_Wh = 3600              # 1 Wh = 3600 J
    J_per_kWh = 3.6e6            # 1 kWh = 3.6e6 J
    J_per_kJ = 1000              # 1 kJ = 1000 J

    # --- To Joules ---
    def eV_to_J(self, eV): return eV * self.J_per_eV
    def cal_to_J(self, cal): return cal * self.J_per_cal
    def kcal_to_J(self, kcal): return kcal * self.J_per_kcal
    def Wh_to_J(self, Wh): return Wh * self.J_per_Wh
    def kWh_to_J(self, kWh): return kWh * self.J_per_kWh
    def kJ_to_J(self, kJ): return kJ * self.J_per_kJ

    # --- From Joules ---
    def J_to_eV(self, J): return J / self.J_per_eV
    def J_to_cal(self, J): return J / self.J_per_cal
    def J_to_kcal(self, J): return J / self.J_per_kcal
    def J_to_Wh(self, J): return J / self.J_per_Wh
    def J_to_kWh(self, J): return J / self.J_per_kWh
    def J_to_kJ(self, J): return J / self.J_per_kJ


class GibbsFreeEnergy(Formula):
    """
    **Gibbs Free Energy**::

        ΔG = ΔH - TΔS

        G : Gibbs free energy (kJ/mol)
        H : Enthalpy          (kJ/mol)
        T : Temperature       (K)
        S : Entropy           (J/(mol·K))
    """
    G: TYPE
    H: TYPE
    T: TYPE
    S: TYPE

    def __init__(self, G=None, H=None, T=None, S=None):
        super().__init__(G, H, T, S)

        G, H, T, S = self.symbols

        self.formula = G - (H - T*S)


class FreeEnergy(Formula):
    """
    **Free Energy**::

        ΔG = -R·T·ln(K_p)

        G   : Gibbs free energy    (kJ/mol)
        R   : Ideal gas constant   (J/(mol·K))
        T   : Temperature          (K)
        K_p : Equilibrium constant
    """
    G: TYPE
    R: TYPE
    T: TYPE
    K_p: TYPE

    def __init__(self, G=None, R=Constants.R_J, T=None, K_p=None):
        super().__init__(G, R, T, K_p)

        G, R, T, K_p = self.symbols

        self.formula = G - (R * T * sp.log(K_p))


class PhotonEnergy(Formula):
    """
    **Photon Energy**::

        E = hv

        E : Photon energy    (J)
        h : Plancks constant (J·s)
        v : Photon frequency (1/s)
    """
    E: TYPE
    h: TYPE
    v: TYPE

    def __init__(self, E=None, h=Constants.h, v=None):
        symbols = sp.symbols('E h v')
        values = [E, h, v]
        self._vars = {name: Variable(sym, val) for name, sym, val in
                      zip(['E', 'h', 'v'], symbols, values)}

        v = self._vars
        self.formula = v['E'] - v['h']*v['v']


class PhotonConstant(Formula):
    """
    **Photon Constant**::

        c = λv

        c : Speed of light (m/s)
        λ : Wave length    (m)
        v : Frequency      (1/s)
    """
    c: TYPE
    l: TYPE
    v: TYPE

    def __init__(self, c=Constants.c, l=None, v=None):
        symbols = sp.symbols('c l v')
        values = [c, l, v]
        self._vars = {name: Variable(sym, val) for name, sym, val in
                      zip(['c', 'l', 'v'], symbols, values)}

        v = self._vars
        self.formula = v['c'] - v['l']*v['v']


class ChangeInBoilingFreezing(Formula):
    """
    **Change in Boiling Freezing**::

        ΔT = i·m·K

        T : Temperature               (K)
        i : Hoff-factor
        m : Molality                  (mol/kg)
        K : Boiling Freezing constant
    """
    T: TYPE
    i: TYPE
    m: TYPE
    K: TYPE

    def __init__(self, T=None, i=1, m=None, K=None):
        symbols = sp.symbols('T i m K')
        values = [T, i, m, K]
        self._vars = {name: Variable(sym, val) for name, sym, val in
                      zip(['T', 'i', 'm', 'K'], symbols, values)}

        v = self._vars
        self.formula = v['T'] - v['i']*v['m']*v['K']


class OsmoticPressure(Formula):
    """
    **Osmotic Pressure**::

        Π = i·M·R·T

        Π : Osmotic Pressure   (atm)
        i : Hoff-factor
        M : Molarity           (mol/L)
        R : Ideal gas constant (L·atm/(mol·K))
        T : Temperature        (K)
    """
    Pi: TYPE
    i: TYPE
    M: TYPE
    R: TYPE
    T: TYPE

    def __init__(self, Pi=None, i=1, M=None, R=Constants.R_L, T=None):
        symbols = sp.symbols('Pi i M R T')
        values = [Pi, i, M, R, T]
        self._vars = {name: Variable(sym, val) for name, sym, val in
                      zip(['Pi', 'i', 'M', 'R', 'T'], symbols, values)}

        v = self._vars
        self.formula = v['Pi'] - v['i']*v['M']*v['R']*v['T']


class Rate(Formula):
    """
    **Rate**::

        rate = [A]^m · [B]^n

        Rate  : Reaction speed (M/s)
        []    : Concentration  (mol/L)
        m | n : Coefficients
    """
    rate: TYPE
    A: TYPE
    m: TYPE
    B: TYPE
    n: TYPE

    def __init__(self, rate=None, A=None, m=None, B=None, n=None):
        symbols = sp.symbols('rate A m B n')
        values = [rate, A, m, B, n]
        self._vars = {name: Variable(sym, val) for name, sym, val in
                      zip(['rate', 'A', 'm', 'B', 'n'], symbols, values)}

        v = self._vars
        self.formula = v['rate'] - v['A']**v['m'] * v['B']**v['n']
        

class FirstOrderConcentration(Formula):
    """
    **First Order Concentration**::

        ln[A] = -kt + ln[A]_0

        []   : Concentration           (mol/L)
        k    : Reaction speed constant (1/s)
        t    : Time                    (s)
        []_0 : Start concentration     (mol/L)
    """
    A: TYPE
    k: TYPE
    t: TYPE
    A_0: TYPE

    def __init__(self, A=None, k=None, t=None, A_0=None):
        symbols = sp.symbols('A k t A_0')
        values = [A, k, t, A_0]
        self._vars = {name: Variable(sym, val) for name, sym, val in
                      zip(['A', 'k', 't', 'A_0'], symbols, values)}

        v = self._vars
        self.formula = sp.log(v['A'].symbol) - (-1 * v['k']*v['t'] + sp.log(v['A_0'].symbol))


class HalfTime(Formula):
    """
    **Half Time**::

        t_1/2 = 0.693 / k
    """
    t_12: TYPE
    k: TYPE

    def __init__(self, t_12=None, k=None):
        symbols = sp.symbols('t_12 k')
        values = [t_12, k]
        self._vars = {name: Variable(sym, val) for name, sym, val in
                      zip(['t_12', 'k'], symbols, values)}

        v = self._vars
        self.formula = v['t_12'] - 0.693 / v['k']


class SecondOrderConcentration(Formula):
    """
    **Second Order Concentration**::

        1/[A] = kt + 1/[A]_0

        []   : Concentration           (mol/L)
        k    : Reaction speed constant (1/(s·M))
        t    : Time                    (s)
        []_0 : Start concentration     (mol/L)
    """
    A: TYPE
    k: TYPE
    t: TYPE
    A_0: TYPE

    def __init__(self, A=None, k=None, t=None, A_0=None):
        symbols = sp.symbols('A k t A_0')
        values = [A, k, t, A_0]
        self._vars = {name: Variable(sym, val) for name, sym, val in
                      zip(['A', 'k', 't', 'A_0'], symbols, values)}

        v = self._vars
        self.formula = 1/v['A'] - (v['k']*v['t'] + 1/v['A_0'])


class ArrheniusDifference(Formula):
    """
    **Arrhenius**::

        ln(k_2/k_1) = - E_a/R · (1/T_2 - 1/T_1)

        k   : Reaction speed constant (1/s)
        E_a : Activation energy       (kJ/mol)
        R   : Ideal gas constant      (J/(mol·K))
        T   : Temperature             (K)
    """
    k_2: TYPE
    k_1: TYPE
    E_a: TYPE
    R: TYPE
    T_1: TYPE
    T_2: TYPE

    def __init__(self, k_2=None, k_1=None, E_a=None, R=Constants.R_J, T_1=None, T_2=None):
        symbols = sp.symbols('k_2 k_1 E_a R T_1 T_2')
        values = [k_2, k_1, E_a, R, T_1, T_2]
        self._vars = {name: Variable(sym, val) for name, sym, val in
                      zip(['k_2', 'k_1', 'E_a', 'R', 'T_1', 'T_2'], symbols, values)}

        v = self._vars
        self.formula = sp.log(v['k_2'].symbol / v['k_1'].symbol) - (-1 * v['E_a']/v['R'] * (1/v['T_2'] - 1/v['T_1']))


class Arrhenius(Formula):
    """
    **Arrhenius**::

        ln(k) = - E_a/R · (1/T) + ln(A)

        k   : Reaction speed constant (1/s)
        E_a : Activation energy       (kJ/mol)
        R   : Ideal gas constant      (J/(mol·K))
        T   : Temperature             (K)
        A   : Frequency factor        (1/s)
    """
    k: TYPE
    E_a: TYPE
    R: TYPE
    T: TYPE
    A: TYPE

    def __init__(self, k=None, E_a=None, R=Constants.R_J, T=None, A=None):
        symbols = sp.symbols('k E_a R T A')
        values = [k, E_a, R, T, A]
        self._vars = {name: Variable(sym, val) for name, sym, val in
                      zip(['k', 'E_a', 'R', 'T', 'A'], symbols, values)}

        v = self._vars
        self.formula = sp.log(v['k'].symbol) - (-1 * v['E_a']/v['R'] * (1/v['T']) + sp.log(v['A'].symbol))


class Equilibrium(Formula):
    """
    **Equilibrium**::

        K = ([C]^c · [D]^d) / ([A]^a · [B]^b)
    """
    K: TYPE
    A: TYPE
    a: TYPE
    B: TYPE
    b: TYPE
    C: TYPE
    c: TYPE
    D: TYPE
    d: TYPE

    def __init__(self, K=None, A=None, a=1, B=None, b=1, C=None, c=1, D=None, d=1):
        symbols = sp.symbols('K A a B b C c D d')
        values = [K, A, a, B, b, C, c, D, d]
        self._vars = {name: Variable(sym, val) for name, sym, val in
                      zip(['K', 'A', 'a', 'B', 'b', 'C', 'c', 'D', 'd'], symbols, values)}

        v = self._vars
        self.formula = v['K'] - (v['C']**v['c'] * v['D']**v['d']) / (v['A']**v['a'] * v['B']**v['b'])


class RadioactiveDecay(Formula):
    """
    **Radioactive Decay**::

        N_t = N_0 · (1/2)^(t/t_1/2)
    """
    N_t: TYPE
    N_0: TYPE
    t: TYPE
    t_12: TYPE

    def __init__(self, N_t=None, N_0=None, t=None, t_12=None):
        symbols = sp.symbols('N_t N_0 t t_12')
        values = [N_t, N_0, t, t_12]
        self._vars = {name: Variable(sym, val) for name, sym, val in
                      zip(['N_t', 'N_0', 't', 't_12'], symbols, values)}

        v = self._vars
        self.formula = v['N_t'] - v['N_0'] * (1/2)**(v['t']/v['t_12'])


class AqueousSolution(Formula):
    pass


class Autoionization(Formula):
    """
    **Autoionization**::

        K_w = [H+][OH-] = 1e-14

        K_w : Equilibrium constant
        []  : Concentration        (M/L)
    """
    K_w: TYPE
    H: TYPE
    OH: TYPE

    def __init__(self, K_w=1e-14, H=None, OH=None):
        super().__init__(K_w, H, OH)

        K_w, H, OH = self.symbols

        self.formula = K_w - H * OH


class PAutoionization(Formula):
    """
    **Autoionization**::

        14 = pH + pOH

        K_w : Equilibrium constant
        []  : Concentration        (M/L)
    """
    pH: TYPE
    pOH: TYPE

    def __init__(self, pH=None, pOH=None):
        super().__init__(pH, pOH)

        pH, pOH = self.symbols

        self.formula = 14 - (pH + pOH)


class PH(Formula):
    """
    **pH**::

        pH = -log([H3O+])
    """
    pH: TYPE
    H3O: TYPE

    def __init__(self, pH=None, H3O=None):
        super().__init__(pH, H3O)

        pH, H3O = self.symbols

        self.formula = pH - (-1 * sp.log(H3O, 10))


class POH(Formula):
    """
    **pH**::

        pOH = -log([OH-])
    """
    pOH: TYPE
    OH: TYPE

    def __init__(self, pOH=None, OH=None):
        super().__init__(pOH, OH)

        pOH, OH = self.symbols

        self.formula = pOH - (-1 * sp.log(OH, 10))


class Buffer(Formula):
    """
    **Buffer solution**::

        pH = pK_a + log([base] / [acid])
    """
    pH: TYPE
    pK_a: TYPE
    base: TYPE
    acid: TYPE

    def __init__(self, pH=None, pK_a=None, base=None, acid=None):
        super().__init__(pH, pK_a, base, acid)

        pH, pK_a, base, acid = self.symbols

        self.formula = pH - (pK_a + sp.log(base/acid, 10))


class CellPotential(Formula):
    """
    **Cell potential**::

        ΔG = -nFE
    """
    G: TYPE
    n: TYPE
    F: TYPE
    E: TYPE

    def __init__(self, G=None, n=None, F=Constants.F, E=None):
        super().__init__(G, n, F, E)

        G, n, F, E = self.symbols

        self.formula = G + n*F*E


class StandardCellPotential(Formula):
    """
    **Standard cell potential**::

        E_cell = E_katode/red - E_anode/ox

        ox: -> e^-
        red: e^- ->
    """
    E_cell: TYPE
    E_katode: TYPE
    E_anode: TYPE

    def __init__(self, E_cell=None, E_katode=None, E_anode=None):
        super().__init__(E_cell, E_katode, E_anode)

        E_cell, E_katode, E_anode = self.symbols

        self.formula = E_cell - (E_katode - E_anode)


class NernstEquation(Formula):
    """
    E = E° - (RT) / (nF) ln(Q)
    """
    E: TYPE
    E_o: TYPE
    R: TYPE
    T: TYPE
    n: TYPE
    F: TYPE
    Q: TYPE

    def __init__(self, E=None, E_o=None, R=Constants.R_J, T=None, n=None, F=Constants.F, Q=None):
        super().__init__(E, E_o, R, T, n, F, Q)

        E, E_o, R, T, n, F, Q = self.symbols

        self.formula = E - (E_o - R*T / (n*F) * sp.log(Q, 10))
