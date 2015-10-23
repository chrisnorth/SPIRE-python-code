from pylab import *
from scipy.optimize import leastsq

iterations = range(1, 31)
time = [454, 163, 167, 496, 199, 214, 528, 240, 541, 273, 268, 303, 312, 321, 336, 348, 357, 382, 384, 403, 406, 427, 453, 457, 467, 479, 499, 503, 528, 548]

iterations = iterations[10:]
time = time[10:]

loglog(iterations, time, 'kx', label='Data')
xlabel('Number of Iterations')
ylabel('Runtime / s')

powerlaw = lambda x, amp, index: amp * (x**index)

logx = log10(iterations)
logy = log10(time)

fitfunc = lambda p, x: p[0] + p[1] * x
errfunc = lambda p, x, y: (y - fitfunc(p, x))
fitfunc2 = lambda p, x: p[0] + 2 * x
errfunc2 = lambda p, x, y: (y - fitfunc2(p, x))

pinit = [1.0, 2.0]
coeffs, _ = leastsq(errfunc, pinit, args=(logx, logy))
coeffs2, _ = leastsq(errfunc2, [1.0], args=(logx, logy))
x_lim = xlim()
y_lim = ylim()

loglog(x_lim, powerlaw(x_lim, 10**coeffs[0], coeffs[1]), 'r-', label='$y\propto x^{{ {:.1f} }}$'.format(coeffs[1]))

# x = logspace(1.5,2.5 , num=200)
# loglog(x, y)
xlim(0,40)
ylim(y_lim)
legend(loc=2)
grid(True, which="both")
savefig('doc/iter-runtime.pdf')
show()

show()
