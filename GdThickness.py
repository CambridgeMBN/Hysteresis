# Area of film

length = 1 # in cm
width = 0.5 # in cm

area = length * width # in cm2

# saturation magnetisation
m_sat = 5.75E-4

# bohr magneton
mu_b =  9.2741E-21 # erg/G, 1 erg/G = 1 emu

# 7.55 bohr magnetom per Gd atom
n_gd = 7.55

# density of Gd = 7.9 g /cm3 or 7.9E3 kg per metre cubed
rho = 7.9E3 # kg / m3

# Atomic weight = 157.25
weight = 157.25

# Unit atomic mass (mass of proton or neutron)
m_unit = 1.67E-27 # kg

# Magnetisation (per unit volume)
m = mu_b * n_gd * (rho / (weight * m_unit)) # in m-3

print 'Number: ', (rho / (weight * m_unit)) # in m-3

# Convert to /cm3
m_cm = m*1E-6
print 'M: ', m_cm
print 'MSat/M: ', m_sat / m_cm
print 'T: ', (m_sat / m_cm) / area

# Calculate thickness
thickness = m_sat / (m_cm * area)

# Thickness is in cm, want in nm 
thickness_nm = thickness * 1E7# *1E9

print 'Thickness: ', thickness_nm