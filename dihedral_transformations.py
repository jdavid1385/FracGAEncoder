# coding: utf-8
from PIL import Image, ImageDraw
from numpy import array, fliplr,flipud,rot90


orig = lambda x: x
refl_y = fliplr
refl_x = flipud
im90 = rot90
im180 = lambda im: rot90(im,2)
im270 = lambda im: rot90(im,3)

transformations = [[orig], [refl_y], [refl_x], [im180], [flipud,im90], [im90], [im270], [flipud,im270]]

def get_dihedrical_transf_all(orig):
    refl_y = fliplr(orig)
    refl_x = flipud(orig)
    im90 = rot90(orig)
    im180 = rot90(im90)
    im45 = flipud(im90)
    im270 = rot90(im180)
    im_min45 = flipud(im270)
    return [orig, refl_y, refl_x, im90, im180, im45, im270, im_min45]

def get_dihedrical_transf(im,transf):
    composeTransf = transformations[transf]
    composeTransf.reverse()
    iterTrans = iter(composeTransf)
    # First have to be a another reference diferent than image
    out = iterTrans.next()(im)
    for f in iterTrans:
        out = f(out)
    composeTransf.reverse()
    return out
