import numpy as np

def separate(x, y, x_cutoff=1.):
    i_start = min(np.argmin(x), np.argmax(x))
    x = np.roll(x, i_start)
    y = np.roll(y, i_start)
    i_te = np.argmax(x)
    x_upper = []
    y_upper = []
    for i in range(1, i_te):
        if x[i] > x[i-1]:
            break
    while i <= i_te and x[i] > x[i-1]:
        x_upper = x_upper + [x[i]]
        y_upper = y_upper + [y[i]]
        i += 1
    x_upper = np.array(x_upper)
    y_upper = np.array(y_upper)
    x_lower = []
    y_lower = []
    for i in range(i_te, x.size):
        if x[i] < x[i-1]:
            break
    while i < x.size and x[i] < x[i-1]:
        x_lower = x_lower + [x[i]]
        y_lower = y_lower + [y[i]]
        i += 1
    x_lower = np.array(x_lower)
    y_lower = np.array(y_lower)
    i = np.argmin(np.abs(x_lower - x_cutoff))
    x_lower = x_lower[i:]
    y_lower = y_lower[i:]    
    i = np.argmin(np.abs(x_upper - x_cutoff))
    x_upper = x_upper[:i]
    y_upper = y_upper[:i]        
    return x_lower, y_lower, x_upper, y_upper
