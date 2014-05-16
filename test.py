# IMPORTS

import signal, os

from numpy import array
from numpy import add

from numpy.random import RandomState
from TransformationCodes import rgb2hex
from heapq import heappush, heappop, nlargest
from heapq import merge as merge_heaps
import scipy.ndimage as get_scaled
from PIL import Image, ImageDraw
from utils import chromosome_init,calculate_mse,hex2rgb,calculate_mse_hex,resize_func,get_phenotype,show_chrom,calc_massic,gen_mask_crossover
from Servers import Generator,Mutator,Fitness_evaluator,Fractal_encoder     
  

# Size of the image
h = 100; w = 100
h_dom = 16; w_dom = 16
N = 10

# GENESIS

if __name__ == '__main__':    
    
    im = Image.new('RGBA', (80, 80), (0, 0, 0, 0)) # Create a blank image

    draw = ImageDraw.Draw(im) # Create a draw object
    draw.rectangle(( 0, 0,  30,  30), fill="blue", outline="yellow")
    draw.rectangle(( 0,50,  30,  80), fill="red", outline="red")
    draw.rectangle((50, 0,  80,  30), fill="green", outline="red")
    Dom = array(im,'ubyte')
    DomBlockSize = (16,16)
    RanBlockSize = (8,8)
    
    N_population=10
    N_workers=2
      
    props_generator = dict( N = N_population, h = 80, w = 80, P_c = 0.65, h_dom = 16, w_dom = 16)
    props_mutator = dict(P_mb=0.005, P_mw = 0.002,  h = 80, w = 80, h_dom = 16, w_dom = 16)
   
    encoder = Fractal_encoder(im, N_population, N_workers, DomBlockSize, RanBlockSize, props_generator, props_mutator)

    Solution = []
    counter = 0
    init_chrom_pool = chromosome_init(N_population, 16, 16, 80, 80)

    for RanBlk in encoder.next_rangBlk():
        
        encoder.update_pool(init_chrom_pool)    
    
        # Counter for number of iterations selecting the same best chromosome
        count_best = 0
        best_individual = ''
        n_iterations = 200

        # 200 iterations at most
        for i in xrange(200):
        
            partition = [(0, 5),(5, 10)]
            mating_pool = encoder.calculate_fitness(partition)

            #S_b = [chrom for meas,chrom in mating_pool]
            
            partition = [(0,3),(3,5)]
            S_w = encoder.generate_new_offsprings(partition)
            encoder.update_pool()
            
            # Stop criteria for the best individual selected up to 16 times
            if best_individual ==  mating_pool[0][1]:
                count_best += 1
                if count_best > 12:
                    n_iterations=i
                    break
            else:
                count_best = 0
                best_individual = mating_pool[0][1]
            
            if i == 6:
                encoder.update_mutation_rate(0.003, 0.002)

            if i == 10:
                encoder.update_mutation_rate(0.001, 0.001)
          
            mutated_population = encoder.apply_mutation(partition)
            encoder.update_pool(mutated_population)
 
        rgb_pq = calc_massic(get_phenotype(best_individual, Dom, DomBlockSize, RanBlockSize ),RanBlk)  
        Solution.append((best_individual, rgb_pq, n_iterations))
        if counter > 10:
            break
        else:
            counter+=1        


    encoder.finish_all()
    print Solution
