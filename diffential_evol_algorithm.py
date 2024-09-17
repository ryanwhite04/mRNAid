import numpy as np
from scipy.optimize import differential_evolution

len = 10

y = [1.3,2.4,3.5,4.1,3.2,4.5,3.7,4.5,6.1,8.2]
# x should be an array
def objective_function(x):
    sum = 0
    for i in range(0,10):
        sum += (x[i]-y[i]) * (x[i]-y[i])

    return sum 
        

# range, not suitable for discrete value
bounds = [(-10,10)]*10

result = differential_evolution(objective_function, bounds)

#print results
print("best solution:", result.x)  
print("min value", result.fun)
print("success:", result.success)
print("exit state:", result.message)
