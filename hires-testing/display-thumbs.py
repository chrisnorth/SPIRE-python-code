from pylab import *
from astropy.io import fits
import aplpy
import glob

#files = glob.glob('/Users/Tom/HIPE/plots/obs_thumbs/above_lowsnr/*PLW.fits')[:16]
# files = glob.glob('/Users/Tom/HIPE/plots/obs_thumbs/below_highsnr/*PLW.fits')[:16]
#
# for f in files:
#     print(f)
#
# fmain = figure(figsize=(12,8))
# for i, f in enumerate(files):
#     print(i, f)
#     fig = aplpy.FITSFigure(f, figure=fmain, subplot=(4,4,i+1))
#     fig.show_colorscale(cmap='gist_heat')
#     fig.axis_labels.hide()
#     fig.tick_labels.hide()
#     fig.add_colorbar()
#     #fig.colorbar.show()
#
#
# tight_layout()
# subplots_adjust(hspace=0.2, wspace=0)
#
# #savefig('doc/above_lowsnr.pdf')
# savefig('doc/below_highsnr.pdf')

# iters = [5,10,15,20,25,30, 35, 40, 50]
# fmain = figure(figsize=(10,8))
# for i, n in enumerate(iters):
#     fig = aplpy.FITSFigure('/Users/Tom/HIPE/plots/iter_test/PLW_%d.fits' % n, figure=fmain, subplot=(3,3,i+1))
#     title(n)
#     fig.show_colorscale(cmap='gist_heat')
#     fig.axis_labels.hide()
#     fig.tick_labels.hide()
#     fig.add_colorbar()
#     fig.recenter(76.85, -6.20, radius=0.12)
#
# tight_layout()
# savefig('doc/iter-thumbs-zoom.pdf')

# files = glob.glob('/Users/Tom/HIPE/plots/obs_thumbs/nearby/*PLW.fits')
# fmain = figure(figsize=(10,8))
# for i, f in enumerate(files):
#     fig = aplpy.FITSFigure(f, figure=fmain, subplot=(3,3,i+1))
#     fig.show_colorscale(cmap='gist_heat')
#     fig.axis_labels.hide()
#     fig.tick_labels.hide()
#     fig.add_colorbar()
#
# tight_layout()
# savefig('doc/nearby-thumbs.pdf')

base_dir = '/Users/Tom/HIPE/plots/radial-beams/'
obsids = [1342249237, 1342227726, 1342210936, 1342216940]
band = 'PLW'

fmain = figure()
for i, obsid in enumerate(obsids):
    hires = base_dir + str(obsid) +'_HIRES_' + band + '.fits'
    radial= base_dir + str(obsid) +'_HIRES_RAIDAL_' + band + '.fits'
    diff = base_dir + str(obsid) +'_DIFF_' + band + '.fits'

    hires_data = fits.open(hires)[1].data
    diff_data = fits.open(diff)[1].data

    fig = aplpy.FITSFigure(hires, figure=fmain, subplot=(4,3,i*3 + 1))
    fig.show_colorscale(cmap='gist_heat')
    fig.axis_labels.hide()
    fig.tick_labels.hide()
    fig.add_colorbar()

    fig = aplpy.FITSFigure(radial, figure=fmain, subplot=(4,3,i*3 + 2))
    fig.show_colorscale(cmap='gist_heat')
    fig.axis_labels.hide()
    fig.tick_labels.hide()
    fig.add_colorbar()

    fig = aplpy.FITSFigure(diff_data / hires_data, figure=fmain, subplot=(4,3,i*3 + 3))
    fig.show_colorscale(cmap='gist_heat')
    fig.axis_labels.hide()
    fig.tick_labels.hide()
    fig.add_colorbar()

# base_dir = '/Users/Tom/HIPE/plots/rotations/'
# obsid = 1342249237
# angles = [30, 60, 90, 120, 150, 180]
# band = 'PLW'
#
# fmain = figure()
# zero_data = fits.open(base_dir + str(obsid) +'_HIRES_ROT_0_' + band + '.fits')[1].data

# for i, angle in enumerate(angles):
#     hires = base_dir + str(obsid) +'_HIRES_ROT_' + str(angle) +  '_' + band + '.fits'
#     diff = base_dir + str(obsid) +'_HIRES_DIFF0_ROT_' + str(angle) +  '_' + band + '.fits'
#     hires_data = fits.open(hires)[1].data
#     diff_data = fits.open(diff)[1].data
#
#     diff_data = diff_data / zero_data
#
#     fig = aplpy.FITSFigure(hires_data, figure=fmain, subplot=(3,4, i*2 + 1))
#     fig.show_colorscale(cmap='gist_heat')
#     fig.axis_labels.hide()
#     fig.tick_labels.hide()
#     fig.add_colorbar()
#     title('Rotated %d°' % angle)
#
#     fig = aplpy.FITSFigure(diff_data, figure=fmain, subplot=(3,4, i*2 + 2))
#     fig.show_colorscale(cmap='gist_heat', vmin=-0.3, vmax=0.3)
#     fig.axis_labels.hide()
#     fig.tick_labels.hide()
#     fig.add_colorbar()
#
# tight_layout()

# base_dir = '/Users/Tom/HIPE/plots/iter_truth/'
# obsid = 1342249237
# band = 'PLW'
# iterations = [0, 1, 5, 10, 20, 30, 50]
#
# fmain = figure()
# fig = aplpy.FITSFigure(base_dir + str(obsid) + '_TRUTH_' + band + '.fits', figure=fmain, subplot=(4, 2, 1))
# fig.show_colorscale(cmap='gist_heat', vmin=0, vmax=2)
# fig.axis_labels.hide()
# fig.tick_labels.hide()
# fig.add_colorbar()
# fig.recenter(76.88, -6.3, radius=0.12)
# title('Truth')
#
# for i, n in enumerate(iterations):
#     if n == 0:
#         fig = aplpy.FITSFigure(base_dir + str(obsid) + '_NOMINAL_' + band + '.fits', figure=fmain, subplot=(4, 2, i+2))
#     else:
#         fig = aplpy.FITSFigure(base_dir + str(obsid) + '_HIRES_' + band + '_' + str(n) + '.fits', figure=fmain, subplot=(4, 2, i+2))
#
#     if n==0:
#         fig.show_colorscale(cmap='gist_heat', vmin=0)
#     else:
#         fig.show_colorscale(cmap='gist_heat')
#     fig.axis_labels.hide()
#     fig.tick_labels.hide()
#     fig.add_colorbar()
#     fig.recenter(76.88, -6.3, radius=0.12)
#
#     if n == 0:
#         title('Nominal')
#     else:
#         title('N=%d' % n)

# base_dir = '/Users/Tom/HIPE/plots/iter_truth/'
# obsid = 1342185538
# band = 'PLW'
# iterations = [0, 1, 5, 10, 20, 30, 50]
#
# fmain = figure()
# fig = aplpy.FITSFigure(base_dir + str(obsid) + '_TRUTH_' + band + '.fits', figure=fmain, subplot=(4, 2, 1))
# fig.show_colorscale(cmap='gist_heat', vmin=0)
# fig.axis_labels.hide()
# fig.tick_labels.hide()
# fig.add_colorbar()
# fig.recenter(148.89, 69.06, radius=0.2)
# title('Truth')
#
# for i, n in enumerate(iterations):
#     if n == 0:
#         fig = aplpy.FITSFigure(base_dir + str(obsid) + '_NOMINAL_' + band + '.fits', figure=fmain, subplot=(4, 2, i+2))
#     else:
#         fig = aplpy.FITSFigure(base_dir + str(obsid) + '_HIRES_' + band + '_' + str(n) + '.fits', figure=fmain, subplot=(4, 2, i+2))
#
#     if n==0:
#         fig.show_colorscale(cmap='gist_heat', vmin=0)
#     else:
#         fig.show_colorscale(cmap='gist_heat')
#     fig.axis_labels.hide()
#     fig.tick_labels.hide()
#     fig.add_colorbar()
#     fig.recenter(148.89, 69.06, radius=0.2)
#
#     if n == 0:
#         title('Nominal')
#     else:
#         title('N=%d' % n)
show()
