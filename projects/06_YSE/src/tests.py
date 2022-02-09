from litho import Litho

myLith = Litho()
myLith.print_vals()
age=1.5 # Myrs
zidx=100
myDepth = myLith.get_depth_sflr(age)
myPressure = myLith.get_obp(myDepth)
myTemp = myLith.get_temperature(age)
print("Age: {} Myr\nDepth: {} m\n".format(age,myDepth))
print("Pressure at z = {} km: {} MPa\n".format(myLith.z[zidx]*1e-3,myPressure[zidx]*1e-6))
print("Temp at z = {} km: {} C\n".format(myLith.z[zidx]*1e-3,myTemp[zidx]))