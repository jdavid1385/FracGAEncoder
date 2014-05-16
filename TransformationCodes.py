# coding: utf-8
from PIL import Image, ImageDraw
from numpy import array
""" The image must be in an nparray instance 
    the height and width in the matrix will
    be in the format (height, width) opossite
    to the representation in PIL (width,height) """
def rgb2hex(img):
    h,w,doncare = img.shape    # The dimensions are inverted in the generated matrix
	# Calculate the rgb code from the hexagesimal representation. This applies only to
	# the internal representation for the PIL Library [0x_A_B_G_R]
    r = img[:,:,0].flatten()
    g = img[:,:,1].flatten()
    b = img[:,:,2].flatten()
    a = img[:,:,3].flatten()
    hexImage = array(map(lambda a_ij, b_ij, g_ij, r_ij: (a_ij<<24) + (b_ij << 16) + (g_ij << 8) + r_ij, a,b,g,r),'uint32')
    return hexImage.reshape((h,w))

