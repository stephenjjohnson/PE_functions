# PE_functions.py

Python functions to perform common petroleum engineering calculations:

`sg2api(sg)`
Convert specific gravity at 60F to API gravity

`api2sg(api)`
Convert API gravity to specific gravity at 60 &deg;F

`darcy_perm(mu, L, Q, A, DeltaP)`
Calculate permeability (k, Darcy) from Darcy's law

`app_visc(k, A, DeltaP, L, Q)`
Calculate apparent viscosity (mu_app, cp) from Darcy's law

`flow_rate(A, DeltaP, k, mu, L)`
Calculate flow rate (Q, cm^3/s) from Darcy's law
   
`shear2rpm(shear)`
Convert shear rate in 1/s to rpm for Brookfield CP-40 viscometer cone

`rpm2shear(rpm)`
Convert rpm to shear rate in 1/s for Brookfield CP-40 viscometer cone

`mobility(k, mu)`
Calculate fluid mobility from permeability and viscosity

`RF(mob_brine, mob_polymer)`
Calculate Resistance Factor from mobility of brine and polymer

`RRF(mob_brine0, mob_brine1)`
Calculate Residual Resistance Factor from brine mobility before and after polymer treatment

`falling_head_perm(a, A, L, t, h_0, h_t, mu=0.9982)`
Calculate permeability coefficient (K) and permeability (k, Darcy) for a falling head

`c2f(degrees_C)`
Convert &deg;C to &deg;F

`f2c(degrees_F)`
Convert &deg;F to degrees &deg;C

`atm2psi(atmosphere)`
Convert atmosphere to pounds per square inch (psi)

`psi2atm(psi)`
Convert pounds per square inch (psi) to atmospheres 

`g2mL(grams, density, porosity)`
Convert grams of sand to volume (mL)

`mL2g(bulkvolume, density, porosity)`
Convert mL of sand to grams

`g2lbm(grams)`
Convert g to lbs

`lbm2g(pounds)`
Convert lbs to g

`retention_lbm_per_acreft(g_per_kg, graindensity = 1, porosity = 0)`
Convert retention from g/kg to lbm per acre foot
   
`retention_g_per_kg(lbm_per_acreft, graindensity = 1, porosity = 0)`
Convert retention from lbm per acre foot to g/kg
   
`hlb_griffin(mass_hydrophilic, mass_total)`
 Calculate hydrophilic-lipophilic balance for non-ionic surfactant from ratio of molar mass of hydrophilic group to molar mass of whole molecule.
See Griffin WC (1954),  "Calculation of HLB Values of Non-Ionic Surfactants," Journal of the Society of Cosmetic Chemists 5:249
   
`hlb_davies(n_hydrophilic, value_hydrophilic, n_lipophilic, value_lipophilic)`
Calculate hydrophilic-lipophilic balance for non-ionic surfactant from number and values of hydrophilic and lipophilic groups.
See Davies JT (1957) "A quantitative kinetic theory of emulsion type, I. Physical chemistry of the emulsifying agent," Gas/Liquid and Liquid/Liquid Interface. In Proceedings of the International Congress of Surface Activity 426-438

`hlb_binary(mass_surfA, hlbA, mass_surfB, hlbB)`
Calculate hydrophilic-lipophilic balance for a binary mixture of surfactants

`shear_tube_v(diameter, velocity)`
Calculate shear rate from diameter of a circular tube and velocity
See https://en.wikipedia.org/wiki/Shear_rate
    
`v_tube_from_shear(diameter, shear_rate)`
Calculate velocity required for a given shear rate in a circular tube.
See https://en.wikipedia.org/wiki/Shear_rate
    
`shear_plates_v(h, velocity)`
Calculate shear rate from distance between two parallel plates and velocity
See https://en.wikipedia.org/wiki/Shear_rate

`v_plates_from_shear(diameter, shear_rate)`
Calculate velocity required for a given shear rate  between two parallel plates
See https://en.wikipedia.org/wiki/Shear_rate

`bond(radius, surface_tension, density1=1, density2=0,  g=9.81)`
Bond number of a droplet of fluid (= ratio of gravitational to surface forces)

`cap_num(visc, velo, ift)`
Capillary number (= ratio of viscous to surface forces)


