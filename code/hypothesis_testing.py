# !/usr/env python3

"""Script to begin exploring statistical measures for project data"""

## IMPORTS ##
from stochastic_nutrient_B import *

### FUNCTIONS ###
#### DISTRIBUTIONS #####
def normal_distribution(x, m=0, s=1):
    return 1/(np.sqrt(2*np.pi*s))*np.exp(-(x-m)**2/2*s)


def standard_error(data, freedom=0):
    """return standard error"""
    return np.sqrt(np.var(data)/(len(data)-freedom))

def ttest(m1, m2, freedom=0, two_way=False):
    """Calculate one way or two way ttest between two data sets"""
    se = standard_error(m1, freedom=freedom)
    if two_way:
        se += standard_error(m2, freedom=freedom)
        return (np.mean(m1)-np.mean(m2))/se
    else:
        return (np.mean(m1)-m2)/se

def correlation_coef(x, y):
    """Calculte pearson's correlation coefficient between two values x and y"""
    sx = np.std(x)
    meanx = np.mean(x)
    meany = np.mean(y)
    sy = np.std(y)

    c= [((x[i]-meanx)/sx)*((y[i]-meany)/sy) for i in range(len(x))]
    c = [i for i in c if i > 0 or i < 0 or i == 0]  #omit missing data


    t = ttest(c,0, freedom=1)

    return { "cor": np.sum(c)/(len(c)-1),
             "t": t,
             "p": normal_distribution(t)}
