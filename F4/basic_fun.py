import numpy as np
import matplotlib.pyplot as plt
import datetime
import globals4simulation as tng
from PIL import Image

def rank2pos(rank):
    rank_row, rank_column = rank
    return tng.r0 + tng.r_step*rank_row, tng.c0 + tng.c_step*rank_column


def laplacian(M):
    z_up = np.roll(M, -1, axis=0)
    z_down = np.roll(M, 1, axis=0)
    z_left = np.roll(M, -1, axis=1)
    z_right = np.roll(M, 1, axis=1)
    return (z_up + z_down + z_left + z_right - 4 * M)


def datetime_stamp():
    now = datetime.datetime.now()
    stamp = now.strftime('%y%m%d_%H%M%S')
    return stamp


def show_field(U, fig = None, ax=None, vmin=None, vmax=None, colorbar = False, log=False):
    if log:
        U_show = np.log10(U)
    else:
        U_show = U
    pos = ax.imshow(U_show,
              # interpolation='bilinear',
              # aspect=nr / nc,
              cmap=plt.cm.coolwarm,  # plt.cm.coolwarm,
              vmin=vmin,
              vmax=vmax)  # cmap=plt.cm.Blues,
    ax.set_axis_off()
    if colorbar:
        fig.colorbar(pos, ax=ax)


def show_plot(x, y, ax=None, vmin=None, vmax=None, log = False):
    if log:
        ax.semilogy(x,y)
    else:
        ax.plot(x,y)
    ax.set_ylim(vmin, vmax)


def aHill(low, high, x, n, th):
    r = (x/th)**n
    return low + (high-low) * r / (1+r)


def rHill(low, high, x, n, th):
    return low + (high-low) / (1 + (x/th) ** n)


def gen_circle(canvas, center, radius, set_value=1):
    from math import ceil, floor
    r_center, c_center = center
    for r in range(ceil(r_center - radius), floor(r_center + radius + 1)):
        for c in range(ceil(c_center - radius), floor(c_center + radius + 1)):  # inclusive
            if (r - r_center) ** 2 + (c - c_center) ** 2 <= radius ** 2:
                canvas[r, c] = set_value


def field2fig(myarray, dst_impath, magnitude_range, cmap=plt.cm.gray, inverse = False):
    vmin, vmax = magnitude_range
    if vmin == None:
        vmin = np.nanmin(myarray)
    if vmax == None:
        vmax = np.nanmax(myarray)
    M = (myarray - vmin) / (vmax - vmin)
    M[M > 1] = 1
    M[M < 0] = 0
    if inverse:
        M = 1-M
    im = Image.fromarray(np.uint8(cmap(M) * 255))  # coolwarm, Reds
    im.save(dst_impath)

    return im

