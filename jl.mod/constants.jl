module ChemPhysConst

export Na,R,uR,Vm,M25,M15,h,c,SoL,hc,k,σ,SKconst,g,grav,echarge,me

global const Na = 6.022e23 #mlc mol^-1 (Avogadro's constant)
global const R = uR = 8.314459848 #J K^-1 mol^-1 (universal gas constant)
global const Vm = 22.414e-3 #m^3 mol^-1 (molar volume)
global const M25 = 1013.25Na/298R #mlc cm^-3 (air number density at 15C)
global const M15 = 1013.25Na/288R #mlc cm^-3 (air number density at 25C)
global const h  = 6.626e-34 #Js (Planck's constant)
global const c = SoL = 299792458. #m/s (speed of light)
global const hc = h⋅c #Jm (product of Planck's constant and speed of light)
global const k = 1.381e-23 #J/K (Boltzmann constant)
global const σ = SKconst = 5.671e-8 #W m^-2 K^-4 (Stefan-Boltzmann constant)
global const g = grav = 9.807 #m s^-2 (standard acceleration of gravity)
global const echarge = 1.602e-19 #C (elementary charge)
global const me = 9.109e-31 #kg (electron mass)

end #module ChemPhysConst
