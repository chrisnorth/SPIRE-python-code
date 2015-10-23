from pylab import *
from scipy.optimize import leastsq

beamhsize =  [50, 60, 70, 80, 90, 100, 125, 150, 200, 250]
runtime = [30, 45, 64, 64, 75, 95,124, 187, 418, 569]

loglog(beamhsize, runtime, 'kx', label='Data')
xlabel('Beam half size / pixels')
ylabel('Runtime / s')

powerlaw = lambda x, amp, index: amp * (x**index)

logx = log10(beamhsize)
logy = log10(runtime)

fitfunc = lambda p, x: p[0] + p[1] * x
errfunc = lambda p, x, y: (y - fitfunc(p, x))
fitfunc2 = lambda p, x: p[0] + 2 * x
errfunc2 = lambda p, x, y: (y - fitfunc2(p, x))

pinit = [1.0, 2.0]
coeffs, _ = leastsq(errfunc, pinit, args=(logx, logy))
coeffs2, _ = leastsq(errfunc2, [1.0], args=(logx, logy))
x_lim = xlim()
y_lim = ylim()

x = logspace(1.2, 2.6)
loglog(x, powerlaw(x, 10**coeffs[0], coeffs[1]), 'r-', label='$y\propto x^{{ {:.1f} }}$'.format(coeffs[1]))
loglog(x, powerlaw(x, 10**coeffs2[0], 2), 'b--', label='$\mathcal{{O}}(n^2)$')
# x = logspace(1.5,2.5 , num=200)
# loglog(x, y)
xlim(30,400)
ylim(y_lim)
legend(loc=2)
grid(True, which="both")
savefig('doc/runtime.pdf')
show()
