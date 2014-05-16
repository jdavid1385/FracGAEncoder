from dihedral_transformations import get_dihedrical_transf

from TransformationCodes import rgb2hex
from heapq import heappush, heappop, nlargest
from heapq import merge as merge_heaps
import scipy.ndimage as get_scaled
from PIL import Image, ImageDraw
from multiprocessing import Array
from multiprocessing import Queue
from utils import *
import ctypes, numpy

from numpy import random

# Exception Classes

class MatingError(Exception):
    def __init__(self):
        self.reason=""
    def set_reason(self,reason):
        self.reason=reason

class MSEError(Exception):
    def __init__(self):
        self.reason=""
    def set_reason(self,reason):
        self.reason=reason

    
class Generator(multiprocessing.Process):
    
    def __init__(self, work_queue, result_queue, event_start, props):

        # Metaproccessed
        multiprocessing.Process.__init__(self)
      
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.State = "Wait"
        
        # Domain specific properties
        
        self.mating_pool = props['S_b']
        self.start_calc = event_start
        self.N = props['N']; self.h = props['h'];  self.w = props['w']
        self.P_c = props['P_c']; self.h_dom = props['h_dom']
        self.w_dom = props['w_dom']
        
        self.mask_crossover = random.RandomState()
       
    def run(self):  
        
        count = 0

        while True:
            
            if self.State == 'Wait':
                self.start_calc.wait()
                self.State = "Calculate"
                
            elif self.State == 'Stop':
                break
                
            while self.State == "Calculate":
                
                # get a part
                req = self.work_queue.get()
                   
                if req == "Stop":
                    self.State = "Stop"
            
                elif req == "Wait":
                    self.State = "Wait"
                    self.result_queue.put("ACK")
                    
                else:   
                    try:
                        result = self.mating(req)     
                    except MatingError, E:
                        result = E.reason
                
                    self.result_queue.put(result)
            
    def mating(self, partition):
        offsprings = []
        
        for i in xrange(partition[0],partition[1]):
            
            mating_pool = self.mating_pool[:].split()
            idx_chosen_parent_1 = self.rank_selection()-1
            idx_chosen_parent_2 = self.rank_selection()-1
            
            try:
                parent_1 = mating_pool[idx_chosen_parent_1]
                parent_2 = mating_pool[idx_chosen_parent_2]
            except:
                error = MatingError()
                error.set_reason(("Error_on_mating_pool_maybe_indexes",idx_chosen_parent_1,idx_chosen_parent_2, self.mating_pool))
                raise error
                
            try:
                if random.random()<=self.P_c:
                    offpring_1, offpring_2 = self.crossover(parent_1, parent_2) 
                else:
                    offpring_1, offpring_2 = parent_1, parent_2
            except:
                error = MatingError()
                error.set_reason(("Error_on_crossover",parent1,parent2))
                raise error

            
            offsprings.append(offpring_1)
            offsprings.append(offpring_2)
            
        return offsprings[0:partition[1]-partition[0]]       
    
    def rank_selection(self):
        random.seed()
        i = 1; chosen = False
        while not chosen:
            random.seed()
            potential = random.randint(1, (self.N/2)) 
            random.seed()
            chosen = not random.randint(1, (self.N/2))>potential
        return potential  
    
    def crossover(self, parent_1, parent_2):
    
        child_1 = ''
        child_2 = ''
    
        mask = self.gen_mask_crossover()

        for idx, mask_bit in enumerate(mask):
        
            if mask_bit == '0':
                child_1 += parent_1[idx]
                child_2 += parent_2[idx]
            else:
                child_1 += parent_2[idx]
                child_2 += parent_1[idx]
    
        childs = []
    
        for chromo in [child_1,child_2]:
            childs.append(self.consistency_check(chromo))
            
        return childs
    
    def update_P_crossover(self, P_c):
        self.P_c = P_c
        
    def gen_mask_crossover(self):
        return '{0:021b}'.format(self.mask_crossover.randint(2**21-1))
    
    def consistency_check(self, chrom):
        # out of bounds check
        w, h = self.w,self.h
        w_dom, h_dom = self.w_dom, self.h_dom
        
        gen_x_dom = int(chrom[0:9],2)
        gen_y_dom = int(chrom[9:18],2)
     
        return '{0:09b}'.format(min(max(gen_x_dom,0),w-w_dom))+'{0:09b}'.format(min(max(gen_y_dom,0),h-h_dom))+chrom[18:21]
 

class Mutator(multiprocessing.Process):
    
    def __init__(self, work_queue, result_queue, event_start, props):
        multiprocessing.Process.__init__(self)
        
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
    
        self.pool = props['pool']
        
        self.start_calc = event_start
        self.P_mb = props['P_mb'];   self.P_mw = props['P_mw']
        self.h = props['h'];  self.w = props['w']
        self.h_dom = props['h_dom']; self.w_dom = props['w_dom']
        
        self.State = "Wait"
       
    def run(self):  
        
        count = 0

        while True:
            
            if self.State == 'Wait':
                self.start_calc.wait()
                self.State = "Calculate"
                
            elif self.State == 'Stop':
                break
                
            while self.State == "Calculate":
                
                # get a part
                req = self.work_queue.get()
                   
                if req == "Stop":
                    self.State = "Stop"
            
                elif req == "Wait":
                    self.State = "Wait"
                    self.result_queue.put("ACK")    

                elif req == "Update_mrate":
					# This is another inner state that holds until
                    # the mutation parameters are updated
                    while not self.work_queue.empty(): continue
                    self.result_queue.put("ACK")                        
                    P_m = self.work_queue.get()
                    self.P_mb = P_m[0]
                    self.P_mw = P_m[1]
                    self.result_queue.put("ACK")
                    self.State = "Wait"

                else:                      
                    result = []
                    try:
                        population = self.pool[:].split()
                        N = len(population)
                        best = population[0:N/2]
                        worst = population[N/2:N]
                        for chroma in best[req[0]:req[1]]:
                            result.append(self.mutate_best(chroma))
                        for chroma in worst[req[0]:req[1]]:
                            result.append(self.mutate_worst(chroma))         
                    except:
                        result.append(["Mutation_Error",(population,N)])
                
                    self.result_queue.put(result)
            
    def mutate_best(self, chrom):
        gen_chrom = chrom
        if random.random() <= self.P_mb:
            gen_x_dom_HB = chrom[0:5]
            gen_y_dom_HB = chrom[9:14]
            gen_chrom =  gen_x_dom_HB +'{0:04b}'.format(random.randint(2**4-1)) + gen_y_dom_HB +'{0:04b}'.format(random.randint(2**4-1))+chrom[18:21]
    
        return self.consistency_check(gen_chrom)

    def mutate_worst(self, chrom):
        gen_chrom = chrom
        if random.random() <= self.P_mw:
            gen_x_dom_LB = chrom[5:9]
            gen_y_dom_LB = chrom[14:18]
            gen_chrom = '{0:05b}'.format(random.randint(2**5-1)) + gen_x_dom_LB + '{0:05b}'.format(random.randint(2**5-1)) + gen_x_dom_LB +'{0:03b}'.format(random.randint(2**3-1))
    
        return self.consistency_check(gen_chrom)
    
    def consistency_check(self, chrom):
        # out of bounds check
        w, h = self.w,self.h
        w_dom, h_dom = self.w_dom, self.h_dom
        
        gen_x_dom = int(chrom[0:9],2)
        gen_y_dom = int(chrom[9:18],2)
     
        return '{0:09b}'.format(min(max(gen_x_dom,0),w-w_dom))+'{0:09b}'.format(min(max(gen_y_dom,0),h-h_dom))+chrom[18:21]


class Fitness_evaluator(multiprocessing.Process):
    
    def __init__(self, work_queue, result_queue, event_start, props):
        multiprocessing.Process.__init__(self)
        
        self.start_calc = event_start
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.calc = False
        self.kill_received = False
        
        self.ranBlk = 0
        
        self.h = props['h']
        self.w = props['w']
        
        self.h_dom = props['h_dom']
        self.w_dom = props['w_dom']
        
        self.RanBlockSize = props['RanBlockSize']
        self.pool = props['pool']
        self.Dom = props['Dom']
        
        
        # Initial state
        self.State = "Wait"
        self.ranBlkIter=self.next_rangBlk()
        
    def run(self): 
        
        count = 0

        while True:
            
            if self.State == 'Wait':
                self.start_calc.wait()
                self.State = "Calculate"
                
            elif self.State == 'Stop':
                break
                
            flag_chRanBlk = False
            while self.State == "Calculate":
                
                # get a part
                req = self.work_queue.get()
                
                if req == "Next_ranBlk":
                    if flag_chRanBlk == False:
                        self.set_rangBlk()
                        self.result_queue.put("ACK")                                                
                        flag_chRanBlk = True
                    else:    
                        # I already done it and I was faster I so have to give it back
                        # This will happen until my lazy fitness calculator' collegues 
                        # wake up and do their job
                        self.work_queue.put(req)
                    
                elif req == "Stop":
                    self.State = "Stop"
            
                elif req == "Wait":
                    self.State = "Wait"
                    self.result_queue.put("ACK")      
                    
                else:     
                    result = self.calculate_Population_fitness(self.pool[:].split()[req[0]:req[1]])
                    self.result_queue.put(result)
            
    def next_rangBlk(self):
        for i in xrange (int(self.w/8)):
            for j in xrange (int(self.h/8)):
                ranBlk = self.Dom[j*self.RanBlockSize[0]:j*self.RanBlockSize[0]+ self.RanBlockSize[0], i*self.RanBlockSize[1]:i*self.RanBlockSize[1]+self.RanBlockSize[1]]
                self.ranBlk = rgb2hex(ranBlk.copy())
                yield None   

    def set_rangBlk(self):
        self.ranBlkIter.next()
        
    def calculate_Population_fitness(self, part):
        # (h,w)
        DomBlockSize = (16,16)
        RanBlockSize = (8,8)
        max_fitness = 1000000000
        ranking = []
        
        for chrom in part:
            # genes
            gen_x_dom = chrom[0:9]
            gen_y_dom = chrom[9:18]
            gen_flip  = chrom[18:21]
    
            # fenotypes
            fen_xdom  = int(gen_x_dom,2)  # 2 for binary representation
            fen_ydom  = int(gen_y_dom,2)  
            fen_flip  = int(gen_flip,2)
            
            try:
                DomBlk = self.Dom[fen_ydom:fen_ydom+DomBlockSize[0] ,fen_xdom:fen_xdom+DomBlockSize[1]]
            except:
                return "DomBlkError"
            
            try:
                DomBlk_hex = rgb2hex(DomBlk.copy())
            except:
                return "rgb2hexError"
            
            try:
                temp = get_scaled.geometric_transform(DomBlk_hex, resize_func, output_shape=RanBlockSize)
            except:
                return "transformError"    
            
            try:
                DomBlk_subsampled = get_dihedrical_transf(temp,fen_flip)
            except:
                return "dihedrTransformError"    
            #p,q = calc_massic(DomBlk_subsampled,rngBlk)
            
            try:
                MSE = self.calculate_mse(DomBlk_subsampled)
            except MSEError,E:
                return E.reason
            try:
                rank = min(1/MSE,max_fitness)
            except ZeroDivisionError:
                rank = max_fitness

            heappush(ranking,(rank,chrom))
            
        return ranking     

    def calculate_mse(self, domBlk_subSmpl):
        Error = MSEError()
        
        try:
            im1flat_T = domBlk_subSmpl.flatten()
            im2flat_T = self.ranBlk.flatten()
    
            im1flat = im1flat_T.copy()
            im2flat = im2flat_T.copy()
            length = float(len(im1flat))
        except:
            Error.set_reason(("ran_Blk",self.ranBlk))
            raise Error
        MSE = 0
        
        
        for i in xrange(0,int(length)):
            pix_1 = im1flat[i]
            pix_2 = im2flat[i]
        
            rgb_1 = hex2rgb(pix_1)
            rgb_2 = hex2rgb(pix_2)
            try:
                diff = array(rgb_1)-array(rgb_2)
                MSE += reduce(add, map(lambda k: k**2, diff))
            except:
                Error.set_reason("Reduce")
                raise Error


        MSE =  MSE/(3*length)
        return MSE


from numpy import array

class Fractal_encoder():
    
    def __init__(self, image, N_population, N_workers, DomBlockSize, RanBlockSize, props_generator, props_mutator):
        self.work_queue = Queue()
        self.result_queue = Queue()

        self.first_time=True
        self.generators = []
        self.fit_calculators = []
        self.mutators = []

        self.Dom = array(image,'ubyte')
        self.h, self.w = self.Dom.shape[0:2]
        
        self.DomBlockSize = DomBlockSize
        self.RanBlockSize = RanBlockSize
        
        self.N = N_population
        

        # Syncronized elements for population, best clan and worst clan
        
        dummy_pool = numpy.array(''.zfill(N_population*21+N_population-1)).tostring()
        self.pool = Array('c', dummy_pool)
               
        dummy_pool_1 = numpy.array(''.zfill((N_population/2)*21+(N_population/2)-1)).tostring()
        self.S_b = Array('c', dummy_pool_1)
        self.S_w = Array('c', dummy_pool_1)
        
        props_generator.setdefault('S_b', self.S_b)
        props_mutator.setdefault('pool', self.pool)
        
        self.start_fit_calc = multiprocessing.Event()
        self.start_generators = multiprocessing.Event()
        self.start_mutators = multiprocessing.Event()
        
        
        for gen in xrange(N_workers):
            self.generators.append(Generator(self.work_queue, self.result_queue, self.start_generators,props_generator))
            self.generators[gen].start()
        
        props_fit_calc = dict( h = self.h, w = self.w, h_dom = 16, w_dom = 16, RanBlockSize = self.RanBlockSize, Dom = self.Dom, pool=self.pool)    
        
        for fit in xrange(N_workers):
            self.fit_calculators.append(Fitness_evaluator(self.work_queue, self.result_queue, self.start_fit_calc,props_fit_calc))
            self.fit_calculators[fit].start()
            
        for mut in xrange(N_workers):
            self.mutators.append(Mutator(self.work_queue, self.result_queue, self.start_mutators, props_mutator))
            self.mutators[mut].start()

        self.curr_ranBlk = []
            
    def update_pool(self,population=None):
        blank = ' '
        offset = 0
        if population==None:
            population = self.S_b[:].split()+self.S_w[:].split()
            
        chrom_len = len(population[0])
        N = len(population)
        for idx, td in enumerate(population):
            if not offset+21 == N*chrom_len+N-1:
                self.pool[offset:offset+22] = td+blank
                offset+=21+1
            else:
                self.pool[offset:offset+21] = td
           
                
    def update_Sw(self,population):
        offset = 0
        blank = ' '
        chrom_len = len(population[0])
        N = len(population)
        for idx, td in enumerate(population):
            if not offset+21 == N*chrom_len+N-1:
                self.S_w[offset:offset+22] = td+blank
                offset+=21+1
            else:
                self.S_w[offset:offset+21] = td

    def update_Sb(self,population):
        offset = 0
        blank = ' '
        chrom_len = len(population[0])
        N = len(population)
        for idx, td in enumerate(population):
            if not offset+21 == N*chrom_len+N-1:
                self.S_b[offset:offset+22] = td+blank
                offset+=21+1
            else:
                self.S_b[offset:offset+21] = td
           
    def next_rangBlk(self):
        
        for i in xrange (int(self.w/8)):
            for j in xrange (int(self.h/8)):
                ranBlk = self.Dom[j*self.RanBlockSize[0]:j*self.RanBlockSize[0]+ self.RanBlockSize[0], i*self.RanBlockSize[1]:i*self.RanBlockSize[1]+self.RanBlockSize[1]]
                self.curr_ranBlk = rgb2hex(ranBlk.copy())
                
                for fit_calculator in self.fit_calculators:
                    self.work_queue.put("Next_ranBlk")

                self.start_fit_calc.set()
                
                ACKS = []
                while len(ACKS) < len(self.fit_calculators):
                    ACKS.append(self.result_queue.get())
                    
                self.start_fit_calc.clear()                    
                    
                yield self.curr_ranBlk
        
    def calculate_fitness(self, partition):
        count = 0
        
        for part in partition:
            self.work_queue.put(part)
            
        self.start_fit_calc.set()
        
        # collect the results off the queue
        results = []
        
        while len(results) < len(partition):
            result = self.result_queue.get()
            results.append(result)
        
        queue = []
        for h in results:
            queue = merge_heaps(queue,h)
        
        mating_pool = nlargest(self.N/2,queue)              # Return only the half upper part,
                                                            # thise will be to the superior clan
                                                            # and are chosen chroms for mating    
        self.start_fit_calc.clear()            
        for calculator in self.fit_calculators:
            self.work_queue.put("Wait")
            
        ACKS = []
        while len(ACKS) < len(self.fit_calculators):
            ACKS.append(self.result_queue.get())            
            
        self.update_Sb([chrom for meas,chrom in mating_pool])
        
        return mating_pool        
    
    def generate_new_offsprings(self, partition):

        for part in partition:
            self.work_queue.put(part)
            
        self.start_generators.set()            
        # collect the results off the queue
        results = []
        while len(results) < len(partition):
            result = self.result_queue.get()
            results.append(result)
          
        offsprings = []
        for resulting in results:
            offsprings+=resulting

        self.start_generators.clear()            
        for generators in self.generators:
            self.work_queue.put("Wait")
        
        
        ACKS = []
        while len(ACKS) < len(self.generators):
            ACKS.append(self.result_queue.get())  
            
        self.update_Sw(offsprings)            
        return offsprings
    
    def apply_mutation(self, partition):
        count = 0
        
        for part in partition:
            self.work_queue.put(part)
            
        self.start_mutators.set()    
        
        # collect the results off the queue
        results = []
        while len(results) < len(partition):
            result = self.result_queue.get()
            results.append(result)
        
        mutated_population = []
        for result in results:
            mutated_population+=result

        self.start_mutators.clear()            
        for mutator in self.mutators:
            self.work_queue.put("Wait")

        ACKS = []
        while len(ACKS) < len(self.mutators):
            ACKS.append(self.result_queue.get())             
            
        return mutated_population

    def update_mutation_rate(self, rate_Pmb, rate_Pmw):

        for mutator in self.mutators:
            self.work_queue.put("Update_mrate")

        self.start_mutators.set()
                
        ACKS = []
        while len(ACKS) < len(self.fit_calculators):
            ACKS.append(self.result_queue.get())
                  
        self.start_mutators.clear()
        
        for mutator in self.mutators:
            self.work_queue.put((rate_Pmb, rate_Pmw))

        ACKS = []
        while len(ACKS) < len(self.generators):
            ACKS.append(self.result_queue.get())    


    def finish_all(self):
        
        for fit_calculator in self.fit_calculators:
            self.work_queue.put("Stop")
        
        self.start_fit_calc.set()  
        
        for mutator in self.mutators:
            self.work_queue.put("Stop")
            
        self.start_mutators.set()                  
        
        for generator in self.generators:
            self.work_queue.put("Stop")

        self.start_generators.set()      
        

        for fit_calculator in self.fit_calculators:
            fit_calculator.join()
            
        for mutator in self.mutators:
            mutator.join()
            
        for generator in self.generators:
            generator.join()
            
            
        for fit_calculator in self.fit_calculators:
            fit_calculator.terminate()
            
        for mutator in self.mutators:
            mutator.terminate()
            
        for generator in self.generators:
            generator.terminate()
            
        self.start_generators.clear()                  
        self.start_mutators.clear()      
        self.start_fit_calc.clear()              


