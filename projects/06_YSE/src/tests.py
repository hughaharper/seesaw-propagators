from litho import Litho

myLith = Litho()
myLith.print_vals()
age=0 # Myrs
myDepth = myLith.get_depth_sflr(age)
print("Age: {} Myr\nDepth: {} m".format(age,myDepth))