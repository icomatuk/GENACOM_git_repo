from IPython import get_ipython
get_ipython().magic('reset -sf')

import numpy as np
import matplotlib.pyplot as plt
import os

plt.close('all')

cwd = os.getcwd()

grad_x_fd = [1.1,
             -7.41471503,
             3.02045547,
             -3.63183972e-07]

grad_x_analytic = [1.0,
                   -7.47252747,
                   3.07933824,
                   2.66453526e-15]

grad_y_fd = [16.1,
             0.60601036,
             -0.20636652,
             -6.52973881e-07]

grad_y_analytic = [16.0,
                   0.46703297,
                   -0.19245864,
                   1.77635684e-15]

time_fd = [0.001639799999999969,
           0.0009251999999999594,
           0.0009217999999999726,
           0.0009243999999999364]

time_analytic = [0.0003273999999999777,
                 0.000135200000000113,
                 0.00014249999999993435,
                 0.0001261000000001289]

iteration = np.arange(1,5)

"---------------------------------"
plt.figure(figsize=(7, 5))
plt.plot(iteration, time_fd, '-o', Linewidth=1.5, color='0.0', mew=4.5,  
         markersize=5, fillstyle='none',label='FD approx')
plt.plot(iteration, time_analytic, '--s', Linewidth=1.5, color='r', mew=4.5,  
         markersize=5, fillstyle='none',label='exact')
plt.xlabel("iterations", fontsize=15, name='Arial')
plt.ylabel("Time required (sec)", fontsize=15, name='Arial')
plt.xticks(name='Arial', fontsize=11)
plt.yticks(name='Arial', fontsize=11)
plt.legend(loc='upper right',prop={'size':10})
plt.grid(linestyle='--', linewidth=0.5)
plt.xlim([0.8, 4.2])
"---------------------------------"
plt.figure(figsize=(7, 5))
plt.plot(iteration, grad_x_fd, '-o', Linewidth=1.5, color='0.0', mew=4.5,  
         markersize=5, fillstyle='none',label='FD approx')
plt.plot(iteration, grad_x_analytic, '--s', Linewidth=1.5, color='r', mew=4.5,  
         markersize=5, fillstyle='none',label='exact')

plt.xlabel("iterations", fontsize=15, name='Arial')
plt.ylabel("Gradient over X", fontsize=15, name='Arial')
plt.xticks(name='Arial', fontsize=11)
plt.yticks(name='Arial', fontsize=11)
plt.legend(loc='upper right',prop={'size':10})
plt.grid(linestyle='--', linewidth=0.5)
plt.xlim([0.8, 4.2])
"---------------------------------"
plt.figure(figsize=(7, 5))
plt.plot(iteration, grad_y_fd, '-o', Linewidth=1.5, color='0.0', mew=4.5,  
         markersize=5, fillstyle='none',label='FD approx')
plt.plot(iteration, grad_y_analytic, '--s', Linewidth=1.5, color='r', mew=4.5,  
         markersize=5, fillstyle='none',label='exact')

plt.xlabel("iterations", fontsize=15, name='Arial')
plt.ylabel("Gradient over Y", fontsize=15, name='Arial')
plt.xticks(name='Arial', fontsize=11)
plt.yticks(name='Arial', fontsize=11)
plt.legend(loc='upper right',prop={'size':10})
plt.grid(linestyle='--', linewidth=0.5)
plt.xlim([0.8, 4.2])






















