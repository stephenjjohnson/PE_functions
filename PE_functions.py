# Stephen Johnson
# Functions for oil and gas calculations

#import math

from __future__ import division

def sg2api(sg):
    """
    Convert specific gravity at 60F to API gravity
    """
    api = (141.5/float(sg))-131.5
    return api

def api2sg(api):
    """
    Convert API gravity to specific gravity at 60F
    """
    sg = 141.5/(float(api)+131.5)
    return sg

def darcy_perm(mu, L, Q, A, DeltaP):
    """
    Permeability (k, Darcy) from Darcy's law
    Takes:
    mu = viscosity, cp
    L = length, cm
    Q = flow rate, cm**3/s
    A = Cross-sectional area, cm**2
    DeltaP = Pressure drop, atm

    Returns:
    k = permeability, Darcy
    """
    k = (float(mu)*L*Q)/(A*DeltaP)
    return k


def app_visc(k, A, DeltaP, L, Q):
    """
    Apparent viscosity (mu_app, cp) from Darcy's law
    Takes:
    k = permeability, Darcy
    A = Cross-sectional area, cm**2
    DeltaP = Pressure drop, atm
    L = length, cm
    Q = flow rate, cm**3/s
    
    Returns:
    mu = viscosity, cp
    """
    mu = (float(k)*A*DeltaP)/(L*Q)
    return mu


def flow_rate(A, DeltaP, k, mu, L):
    """
    Flow rate (Q, cm**3/s) from Darcy's law
    
    Takes:
    A = Cross-sectional area, cm**2
    DeltaP = Pressure drop, atm
    k = permeability, Darcy
    mu = viscosity, cp
    L = length, cm
    
    Returns:
    Q = flow rate, cm**3/s
    """
    Q = (float(k)*A*DeltaP)/(mu*L)
    return Q

def shear2rpm(shear):
    """
    Convert shear rate in 1/s to rpm for Brookfield CP-40 viscometer cone
    """
    rpm = float(shear)/7.5
    return rpm

def rpm2shear(rpm):
    """
    Convert rpm to shear rate in 1/s for Brookfield CP-40 viscometer cone
    """
    shear = float(rpm) * 7.5
    return shear

def mobility(k, mu):
    """
    Fluid mobility from permeability and viscosity
    """
    mobility = float(k)/mu
    return mobility

def RF(mob_brine, mob_polymer):
    """
    Resistance Factor
    RF = (k/mu_brine)/(k/mu_polymer)
    """
    RF = float(mob_brine) / mob_polymer
    return RF

def RRF(mob_brine0, mob_brine1):
    """
    Residual Resistance Factor
    RRF = (k/mu_brine)_before-polymer/(k/mu_polymer)_after-polymer
    """
    RRF = float(mob_brine0) / mob_brine1
    return RRF

def falling_head_perm(a, A, L, t, h_0, h_t, mu=0.9982):
    """
    Permeability coefficient (K) and permeability (k, Darcy) for a falling head
    
    Takes:
    a = Cross-sectional area of tubing above filter, cm**2
    A = Cross-sectional area of filter, cm**2
    L = length of filter, cm
    t = time for water level to fall from h_0 to h_t, s
    h_0 = initial head of fluid, cm
    h_t = head of fluid at time = t, cm
    mu = fluid viscosity, cp (defaults to water at 20 degC = 0.9982 g/cm**3)

    Uses:
    Ratio of cm.H2O / atm = 1033.22

    Returns:
    K = permeability coefficient
    k = permeability, Darcy
    """
    K = (float(a)/A) * (L/t) * math.log((h_0/h_t), math.e)
    k = K * 1032.22 * mu
    return (K, k)

def c2f(degrees_C):
    """
    Convert degrees Celsius to degrees Fahrenheit
    """
    degrees_F = ((float(degrees_C) * 9)/5) + 32.0
    return degrees_F

def f2c(degrees_F):
    """
    Convert degrees Fahrenheit to degrees Celsius
    """
    degrees_C = ((float(degrees_F) - 32.0) * 5)/9
    return degrees_C

def atm2psi(atmosphere):
    """
    Convert atmosphere to pounds per square inch (psi)
    """
    psi = 14.6959487755134 * float(atmosphere)
    return psi

def psi2atm(psi):
    """
    Convert pounds per square inch (psi) to atmospheres 
    """
    atmosphere = 0.0680459639 * float(psi)
    return atmosphere

def g2mL(grams, density, porosity):
    """
    Convert grams of sand to volume (mL)
    """
    grainvolume = float(grams) / float(density)
    bulkvolume = grainvolume * (1 - porosity)
    return bulkvolume

def mL2g(bulkvolume, density, porosity):
    """
    Convert mL of sand to grams
    """
    grainvolume = float(bulkvolume) / (1 - porosity)
    grams = float(grainvolume) * float(density)
    return grams

def g2lbm(grams):
    """
    Convert g to lbs
    """
    pounds = 0.002205 * float(grams)
    return pounds

def lbm2g(pounds):
    """
    Convert lbs to g
    """
    grams = 453.6 * float(pounds)
    return grams

def retention_lbm_per_acreft(g_per_kg, graindensity = 1, porosity = 0):
    """
    Convert retention from g/kg to lbm per acre foot
    """
    volume_correction = float(graindensity) / (1-float(porosity))
    g_per_L = g_per_kg  * volume_correction
    lbm_per_acreft = 2831.684659 * float(g_per_L)
    return lbm_per_acreft

def retention_g_per_kg(lbm_per_acreft, graindensity = 1, porosity = 0):
    """
    Convert retention from lbm per acre foot to g/kg
    """
    volume_correction = (1-float(porosity)) / float(graindensity)
    g_per_L = float(lbm_per_acreft) / 2831.684659
    g_per_kg = g_per_L  * volume_correction
    return g_per_kg
    
def hlb_griffin(mass_hydrophilic, mass_total):
    """
    Calculate hydrophilic-lipophilic balance for non-ionic surfactant from ratio of molar mass of hydrophilic group to molar mass of whole molecule.
    Griffin WC: "Calculation of HLB Values of Non-Ionic Surfactants," Journal of the Society of Cosmetic Chemists 5 (1954): 249
    """
    hlb = 20 *(float(mass_hydrophilic)/float(mass_total))
    return hlb

def hlb_davies(n_hydrophilic, value_hydrophilic, n_lipophilic, value_lipophilic):
    """
    Calculate hydrophilic-lipophilic balance for non-ionic surfactant from number and values of hydrophilic and lipophilic groups.
    Davies JT: "A quantitative kinetic theory of emulsion type, I. Physical chemistry of the emulsifying agent," Gas/Liquid and Liquid/Liquid Interface. Proceedings of the International Congress of Surface Activity (1957): 426-438
    """
    hlb = 7 + (float(n_hydrophilic) * float(value_hydrophilic)) - (float(n_lipophilic) * float(value_lipophilic))
    return hlb

def hlb_binary(mass_surfA, hlbA, mass_surfB, hlbB):
    """
    Calculate hydrophilic-lipophilic balance for a binary mixture of surfactants
    """
    total_mass = float(mass_surfA) + float(mass_surfB)
    hlb = ((float(mass_surfA)/total_mass) * float(hlbA)) + ((float(mass_surfB)/total_mass) * float(hlbB))
    return hlb

def shear_tube_v(diameter, velocity):
    """
    Calculate shear rate from diameter of a circular tube and velocity
    Diameter, cm
    Velocity, cm/s
    -> Shear rate, 1/s
    https://en.wikipedia.org/wiki/Shear_rate
    """
    shear_rate = ((8*float(velocity))/diameter)
    return shear_rate

def v_tube_from_shear(diameter, shear_rate):
    """
    Calculate velocity required for a given shear rate in a circular tube
    Diameter, cm
    Shear rate, 1/s
    -> Velocity, cm/s
    https://en.wikipedia.org/wiki/Shear_rate
    """
    velocity = ((float(diameter)*shear_rate)/8)
    
def shear_plates_v(h, velocity):
    """
    Calculate shear rate from distance between two parallel plates and velocity
    h = distance between parallel plates, cm
    Velocity, cm/s
    -> Shear rate, 1/s
    https://en.wikipedia.org/wiki/Shear_rate
    """
    shear_rate = float(velocity)/h
    return shear_rate

def v_plates_from_shear(diameter, shear_rate):
    """
    Calculate velocity required for a given shear rate between two parallel plates
    h = distance between, cm
    Shear rate, 1/s
    -> Velocity, cm/s
    https://en.wikipedia.org/wiki/Shear_rate
    """
    velocity = float(h)*shear_rate
    return velocity

def bond(radius, surface_tension, density1=1, density2=0,  g=9.81):
    """
    Bond number of a droplet of fluid (= ratio of gravitational to surface forces)
    Takes:
    density1 = Density of drop phase, g/cm**3
    density2 = Density of surrounding phase, g/cm**3
    radius = Radius of drop, mm
    surface_tension = Surface or interfacial tension, mJ/m**2
    g = Acceleration due to gravity, 9.81 m/s**2
    
    Returns:
    bond = Bond Number (dimensionless)
    """
    delta_rho = float(density1) - float(density2)
    radius_squared = float(radius) * float(radius)
    bond = (delta_rho*g*radius_squared)/surface_tension
    return bond

def cap_num(visc, velo, ift):
    """
    Capillary number (= ratio of viscous to surface forces)
    Takes:
    visc = Viscosity of mobile phase, Poise
    velo = Characteristic velocity of mobile phase, cm**2/s
    ift = Interfacial tension, dyn/cm
    
    Returns:
    capnum = Capillary Number (dimensionless)
    """
    capnum = (float(visc) * velo) / ift
    return capnum
