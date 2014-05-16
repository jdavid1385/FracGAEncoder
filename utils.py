from heapq import heappush, heappop
from numpy import add
from numpy.random import RandomState
from PIL import Image, ImageDraw
import ctypes
from numpy import array
from TransformationCodes import rgb2hex
import scipy.ndimage as get_scaled
from dihedral_transformations import get_dihedrical_transf

mask_crossover = RandomState()

def chromosome_init(N, h_dom=16, w_dom=16, h=100, w=100):
    x_dom = RandomState() #[0-(w-16-1)]
    y_dom = RandomState() #[0-(h-16-1)]
    flip =  RandomState()  #[0-7]

    chrom_pool = (ctypes.c_char_p * N)()
       
    for i in xrange(N):
        chromosome = ctypes.create_string_buffer(2*9+3)          # 9 bits x_dom, 9 bits y_dom and 3 bits flip
        chromosome.value = '{0:09b}'.format(x_dom.randint(0,w-w_dom))+'{0:09b}'.format(y_dom.randint(0,h-h_dom))+'{0:03b}'.format(flip.randint(0,8))
        chrom_pool[i] = ctypes.string_at(chromosome,21)
    return chrom_pool


def calculate_mse(domBlk_subSmpl,rngBlk):
    lenght = float(len(rngBlk.flatten()))
    for i in xrange(0,30):
        R = rngBlk-domBlk_subSmpl
        diff = R[:,:,0:3].flatten()
        MSE = reduce(add, map(lambda k: k**2, diff))/(3*lenght)

def hex2rgb(hexcolor):
    b = ( hexcolor >> 16 ) & 0xFF
    g = ( hexcolor >> 8  ) & 0xFF
    r = hexcolor & 0xFF
    return [r,g,b]
   
def calculate_mse_hex(domBlk_subSmpl,rngBlk):
    im1flat_T = domBlk_subSmpl.flatten()
    im2flat_T = rngBlk.flatten()
    
    im1flat = im1flat_T.copy()
    im2flat = im2flat_T.copy()
    length = float(len(im1flat))
    
    MSE = 0
    for i in xrange(0,length):
        pix_1 = im1flat[i]
        pix_2 = im2flat[i]
        
        rgb_1 = hex2rgb(pix_1)
        rgb_2 = hex2rgb(pix_2)
        
        diff = array(rgb_1)-array(rgb_2)
        MSE += reduce(add, map(lambda k: k**2, diff))

    MSE =  MSE/(3*length)
    return MSE

def resize_func(output_coords):
    return int(float(output_coords[0])*2), int(float(output_coords[1])*2)

def get_phenotype(chrom, Dom,DomBlockSize, RanBlockSize ):
    gen_x_dom = chrom[0:9]
    gen_y_dom = chrom[9:18]
    gen_flip  = chrom[18:21]
    
    # fenotypes
    fen_xdom  = int(gen_x_dom,2)  # 2 for binary representation
    fen_ydom  = int(gen_y_dom,2)  
    fen_flip  = int(gen_flip,2)
    
    DomBlk = Dom[fen_ydom:fen_ydom+DomBlockSize[0] ,fen_xdom:fen_xdom+DomBlockSize[1]]
    DomBlk_hex = rgb2hex(DomBlk.copy())
    temp = get_scaled.geometric_transform(DomBlk_hex, resize_func, output_shape=RanBlockSize)
    DomBlk_subsampled = get_dihedrical_transf(temp,fen_flip)
    
    return DomBlk_subsampled

def show_chrom(chrom, Dom, DomBlockSize ):

    RanBlockSize = (8,8)
     
    gen_x_dom = chrom[0:9]
    gen_y_dom = chrom[9:18]
    gen_flip  = chrom[18:21]

    # fenotypes
    fen_xdom  = int(gen_x_dom,2)  # 2 for binary representation
    fen_ydom  = int(gen_y_dom,2)  
    fen_flip  = int(gen_flip,2)

    DomBlk = Dom[fen_ydom:fen_ydom+DomBlockSize[0] ,fen_xdom:fen_xdom+DomBlockSize[1]]
    DomBlk_hex = rgb2hex(DomBlk.copy())

    temp = get_scaled.geometric_transform(DomBlk_hex, resize_func, output_shape=RanBlockSize)
    DomBlk_subsampled = get_dihedrical_transf(temp,fen_flip)
    
    DomBlk_subsampled = DomBlk_subsampled.copy()
    pilImage = Image.frombuffer('RGBA',DomBlk_subsampled.shape,DomBlk_subsampled,'raw','RGBA',0,1) #for rgba
    #pilImage = Image.frombuffer('RGBA',DomBlk_hex.shape,DomBlk_hex,'raw','RGBA',0,1) #for rgba

    imshow(pilImage.transpose(1))

def calc_massic(domBlk_subSmpl,rngBlk):
    
    u_Temp = domBlk_subSmpl.flatten()
    v_Temp = rngBlk.flatten()
    
    u = u_Temp.copy()
    v = v_Temp.copy()
    
    N = float(len(v))
    
    r_u = []
    g_u = []
    b_u = []
    
    r_v = []
    g_v = []
    b_v = []
    
    for i in xrange(0,N):
        
        pix_u = u[i]
        pix_v = v[i]
        
        [pix_ru, pix_gu, pix_bu] = hex2rgb(pix_u)
        [pix_rv, pix_gv, pix_bv] = hex2rgb(pix_v)
        
        r_u.append(pix_ru)
        g_u.append(pix_gu)
        b_u.append(pix_bu)
        
        r_v.append(pix_rv)
        g_v.append(pix_gv)
        b_v.append(pix_bv)
    
    cum_sum_r = 0; cum_sum_g = 0;  cum_sum_b = 0
    sum_u_r = 0; sum_u_g = 0; sum_u_b = 0
    sum_v_r = 0; sum_v_g = 0; sum_v_b = 0
    sqr_sum_u_r = 0; sqr_sum_u_g = 0; sqr_sum_u_b = 0
    
    
    for (pix_ru,pix_rv),(pix_gu,pix_gv),(pix_bu,pix_bv) in zip(zip(r_u,r_v),zip(g_u,g_v),zip(b_u,b_v)):
        cum_sum_r += pix_ru*pix_rv
        cum_sum_g += pix_gu*pix_gv
        cum_sum_b += pix_bu*pix_bv
        
        sum_u_r += pix_ru
        sum_u_g += pix_gu
        sum_u_b += pix_bu
        
        sum_v_r += pix_rv
        sum_v_g += pix_gv
        sum_v_b += pix_bv
        
        sqr_sum_u_r += pix_ru**2
        sqr_sum_u_g += pix_gu**2
        sqr_sum_u_b += pix_bu**2
    
    try:
        den = N*sqr_sum_u_r-sum_u_r**2
        if not den == 0:
            p_r = round((N*cum_sum_r - sum_u_r*sum_v_r)/ den)
            q_r = round((sum_v_r-p_r*sum_u_r)/N)
        else:
            p_r = 1
            q_r = 0
    except:
        p_r = 1
        q_r = 0
     
    try:
        den = N*sqr_sum_u_g-sum_u_g**2
        if not den == 0: 
            p_g = round((N*cum_sum_g - sum_u_g*sum_v_g)/den)
            q_g = round((sum_v_g-p_g*sum_u_g)/N)
        else:
            p_g = 1
            q_g = 0
    except:
        p_g = 1
        q_g = 0
    
    try:
        den = N*sqr_sum_u_b-sum_u_b**2
        if not den == 0:
            p_b = round((N*cum_sum_b - sum_u_b*sum_v_b)/ den )         
            q_b = round((sum_v_b-p_b*sum_u_b)/N)    
        else:
            p_b = 1
            q_b = 0
    except:        
        p_b = 1
        q_b = 0            
    
    
    return ((p_r, q_r), (p_g, q_g), (p_b, q_b))

def gen_mask_crossover():
    global mask_crossover
    return '{0:021b}'.format(mask_crossover.randint(2**21-1))
    
