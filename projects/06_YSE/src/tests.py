from litho import Litho

myLith = Litho()
myLith.print_vals()
age=1.5 # Myrs
zidx=250
myDepth = myLith.get_depth_sflr(age)
myPressure = myLith.get_obp(myDepth)
myTemp = myLith.get_temperature(age)
myDuct = myLith.get_ductile(myTemp)
myByerC = myLith.get_byerlee(myPressure,'c')
myByerT = myLith.get_byerlee(myPressure,'t')
myByerS = myLith.get_byerlee(myPressure,'s')

print("Age: {} Myr\nDepth: {} m".format(age,myDepth))
print("Pressure at z = {} km: {} MPa".format(myLith.z[zidx]*1e-3,myPressure[zidx]*1e-6))
print("Temp at z = {} km: {} C".format(myLith.z[zidx]*1e-3,myTemp[zidx]))
print("Ductile Str at z = {} km: {} MPa".format(myLith.z[zidx]*1e-3,myDuct[zidx]*1e-6))
print("Strength in compression = {} MPa".format(myByerC[zidx]*1e-6))
print("Strength in tension = {} MPa".format(myByerT[zidx]*1e-6))
print("Strength in shear = {} MPa".format(myByerS[zidx]*1e-6))