#!/usr/bin/env python
import scipy,random,scipy.linalg,scipy.optimize,scipy.sparse,datetime,copy,pickle,math,pickle,numpy

def f_opt(x,medians,C_post):
    #print 'in f'

    inversee    = numpy.linalg.inv(C_post)
    second_term = numpy.dot(inversee,(x-medians))
   
    return numpy.dot((x-medians).transpose(),second_term)
    
#def fmin_differential_evolution(f,x0,population_size=100,generations=20000,bounds=None,variable_is_logarithmic=None,crossover_factor=0.2,disp=0):
def fmin_gen(f,x0,population_size=100,survivors=20,generations=20000,bounds=None,variable_is_logarithmic=None,intruders=0,use_pp=True,convenience_class=None, disp=0):
    import struct

    f = open('medians.txt','r')
    medians_no = []
    content = f.read()
    for elem in content.split(','):
        medians_no.append(float(elem))
    f.close()
    
    g = open('cpost.txt','r')
    content2 = g.read()
    C_post_no = []
    for line in content2.split('\n'):
        if line != '':
            linecontent = line.split(',')
            single_line = []
            for elem in linecontent:
                single_line.append(float(elem))
            C_post_no.append(single_line)

    medians = numpy.array(medians_no)
    C_post  = numpy.array(C_post_no)

    def local_optimize(indiv,convenience_class=None):
        better_indiv = indiv
        if convenience_class:
            better_indiv = scipy.optimize.fmin_bfgs(convenience_class.f,better_indiv,disp=0,maxiter=20)
            fval = convenience_class.f_opt(better_indiv,medians,C_post)
        else:
            #better_indiv = scipy.optimize.fmin_bfgs(f,better_indiv,disp=0,maxiter=20)
            fval = f_opt(better_indiv,medians,C_post)
        return [better_indiv,fval]
    def floatToBits(value):
        if type(value) == type(1.0):
            return (str(struct.unpack('Q', struct.pack('d', value))[0])).rjust(20,"0")
        else:
            return ",".join([floatToBits(x) for x in value.tolist()])
    def bitsToFloat(bits):
        if "," not in bits:
            return struct.unpack('d', struct.pack('Q', long(bits)))[0]
        else:
            return scipy.array([bitsToFloat(x) for x in bits.split(",")])

    def new_individual():
        x = []
        for i in range(indiv_size):
            if variable_is_logarithmic[i]:
                logmin = scipy.log(bounds[i][0])
                logmax = scipy.log(bounds[i][1])
                x.append(scipy.exp(scipy.rand()*(logmax-logmin)+logmin))
            else:
                x.append(scipy.rand()*(bounds[i][1]-bounds[i][0])+bounds[i][0])
        return scipy.array(x)

    def bool_mate(mother,father):
        ms = floatToBits(mother)
        fs = floatToBits(father)
        l = random.sample(range(len(ms)),3)
        l.sort()
        [i1,i2,i3] = l
        cs = ms[:i1]+fs[i1:i2]+ms[i2:i3]+fs[i3:]
        child = bitsToFloat(cs)
        if child.size != mother.size or None in child or scipy.inf in child or -scipy.inf in child: raise ValueError()
        return child

    def mate(mflist):
        #return bool_mate(mflist[0],mflist[1])
        return mflist[0]+mflist[1]-mflist[2]

    def mutate(indiv):
        if not False:
            bi = floatToBits(indiv)
            number = int(math.ceil(float(len(bi))/100.))
            change_indices = random.sample([x for x in range(3,len(bi))],number)
            for ci in change_indices:
                new = bi[:ci]+str(random.choice(range(10)))+bi[(ci+1):]
                try:
                    bitsToFloat(new)
                    bi = new
                except struct.error:
                    # dont accept change
                    pass
            return scipy.absolute(bitsToFloat(bi))
        else:
            return scipy.exp(scipy.log(indiv)+wolf["mutation_factor"]*scipy.array([random.normalvariate(0,1) for x in range(indiv_size)]))

    def bound(vector):

        global correct_vector

        if len(vector) == indiv_size:
            correct_vector = vector

        if len(vector) != indiv_size:
            print 'im doing it right now!'
            vector = correct_vector
        
        for i in range(indiv_size):
            if vector[i] < bounds[i][0]:
                vector[i] = bounds[i][0]
                correct_vector = vector
            elif vector[i] > bounds[i][1]:
                vector[i] = bounds[i][1]
                correct_vector = vector
        return vector

    if bounds is None:
        bounds = [[1e-4, 1e4]]*len(x0)
    if len(bounds) != len(x0):
        raise Exception('Length of x0 and length of bounds do not fit!')

    if variable_is_logarithmic is None:
        variable_is_logarithmic = [[True]]*len(x0)
    if len(variable_is_logarithmic) != len(x0):
        raise Exception('Length of variable_is_logarithmic and x0 do not fit!')

    if use_pp:
        import pp
        # get pp servers
        servers = []
        fi = file("servers","r")
        for l in fi:
            if l.startswith("#"): continue
            servers.append(l.split(",")[1])
        fi.close()
        # read secret file
        try:
            fi = open("secret","r")
            secret = fi.readlines()[0]
            fi.close()
        except:
            import sys
            print "Please create a file calles \"secret\" with one long string in the first line."
            sys.exit(1)
        if disp==1:
            print "starting servers",servers
        job_server = pp.Server(ppservers=tuple(servers),secret=secret)
        print job_server.get_active_nodes()

    population = [x0]
    indiv_size = x0.size
    quality = []
    for i in range(population_size-len(population)):
        population.append(new_individual())
    best_indiv = copy.deepcopy(x0)
    try:
        for i in range(generations):
            # rate individuals
            local_optimization_results = []
            pre_computed_qualities = len(quality)
            for j in range(pre_computed_qualities,population_size):
                if use_pp:
                    if convenience_class:
                        local_optimization_results.append( job_server.submit(local_optimize,(population[j],convenience_class),modules=("scipy","scipy.optimize","scipy.linalg","bounded"),globals=globals() ) )
                    else:
                        local_optimization_results.append( job_server.submit(local_optimize,(population[j],),modules=("scipy","scipy.optimize","scipy.linalg","bounded"),globals=globals() ) )
                else:
                    local_optimization_results.append( local_optimize(population[j]) )
            if use_pp:
                for j in range(len(local_optimization_results)):
                    local_optimization_results[j] = local_optimization_results[j]()
            for j in range(len(local_optimization_results)):
                [better_indiv,fval] = local_optimization_results[j]
                population[pre_computed_qualities+j] = better_indiv
                quality.append( fval )
            # replace None results
            for j in range(len(population)):
                while quality[j] == None or scipy.isnan(quality[j]):
                    population[j] = new_individual()
                    quality[j] = f_opt(population[j],medians,C_post)

            # sort
            sorted_quality = zip(quality,population)
            sorted_quality.sort(key=lambda x:x[0])
            new_population = []
            new_quality = []
            for j in range(survivors):
                f_val,indiv = sorted_quality[j]
                new_population.append(indiv)
                new_quality.append(f_val)
            population = new_population
            quality = new_quality
            # intrude
            for j in range(intruders):
                population.append(new_individual())
            # mate
            current_size = len(population)
            for j in range(current_size,population_size):
                population.append(bound(mutate(mate(random.sample(population[:current_size],3)))))
            best_indiv = copy.deepcopy(population[0])
            if disp==1:
                print "generation",str(i+1).rjust(7),"     f =",quality[0]
                print best_indiv
    except KeyboardInterrupt:
        if disp==1:
            print "Stopping computation"
    if use_pp:
        job_server.print_stats()
    if disp==1:
        print "generation goodbye      f =",quality[0]

    return best_indiv


def fmin_brute(f,ranges,logarithmic=True):
    if type(ranges) == numpy.ndarray:
        newranges = []
        numpoints = 2
        for i in range(len(ranges)):
            newranges.append( [ranges[i]/2.0,ranges[i]*2.0,numpoints] )
        ranges = newranges
    y_to_x = {}
    def range_to_list(tuple):
        if not logarithmic:
            return [((tuple[1]-tuple[0])/tuple[2]*y+tuple[0]) for y in range(tuple[2]+1)]
        else:
            return [math.exp((math.log(tuple[1])-math.log(tuple[0]))/tuple[2]*y+math.log(tuple[0])) for y in range(tuple[2]+1)]
    def fix_part(f,ranges,fixed_x):
        if len(ranges) == 0:
            y = f(scipy.array(fixed_x))
            y_to_x[y] = scipy.array(fixed_x)
        else:
            for this_x in range_to_list(ranges[0]):
                fix_part(f,ranges[1:],fixed_x+[this_x])
    fix_part(f,ranges,[])
    allys = y_to_x.keys()
    allys.sort()
    return y_to_x[allys[0]]
        

def fmin_differential_evolution(f,x0,population_size=100,generations=20000,bounds=None,variable_is_logarithmic=None,crossover_factor=0.2,disp=0):
    import struct

    p = open('medians.txt','r')
    medians_no = []
    content = p.read()
    for elem in content.split(','):
        medians_no.append(float(elem))
    p.close()
    
    g = open('cpost.txt','r')
    content2 = g.read()
    C_post_no = []
    for line in content2.split('\n'):
        if line != '':
            linecontent = line.split(',')
            single_line = []
            for elem in linecontent:
                single_line.append(float(elem))
            C_post_no.append(single_line)

    medians = numpy.array(medians_no)
    C_post  = numpy.array(C_post_no)
    
    def floatToBits(value):
        if type(value) == type(1.0):
            return (str(struct.unpack('Q', struct.pack('d', value))[0])).rjust(20,"0")
        else:
            return ",".join([floatToBits(x) for x in value.tolist()])
    def bitsToFloat(bits):
        if "," not in bits:
            return struct.unpack('d', struct.pack('Q', long(bits)))[0]
        else:
            return scipy.array([bitsToFloat(x) for x in bits.split(",")])
    def new_individual():
        x = []
        for i in range(indiv_size):
            if variable_is_logarithmic[i]:
                logmin = scipy.log(bounds[i][0])
                logmax = scipy.log(bounds[i][1])
                x.append(scipy.exp(scipy.rand()*(logmax-logmin)+logmin))
            else:
                x.append(scipy.rand()*(bounds[i][1]-bounds[i][0])+bounds[i][0])
        return scipy.array(x)
    def crossover(orig,crossing_vector):
        for i in range(len(orig.tolist())):
            if random.random()<crossover_factor:
                orig[i] = crossing_vector[i]
    def bound(vector):
        for i in range(indiv_size):
            if vector[i] < bounds[i][0]:
                vector[i] = bounds[i][0]
            elif vector[i] > bounds[i][1]:
                vector[i] = bounds[i][1]

    if bounds is None:
        bounds = [[1e-4, 1e4]]*len(x0)
    if len(bounds) != len(x0):
        raise Exception('Length of x0 and length of bounds do not fit!')
    if variable_is_logarithmic is None:
        variable_is_logarithmic = [[True]]*len(x0)
    if len(variable_is_logarithmic) != len(x0):
        raise Exception('Length of variable_is_logarithmic and x0 do not fit!')

    population = [x0]
    indiv_size = x0.size
    quality = [f(x0)]
    for i in range(population_size-len(population)):
        population.append(new_individual())
        quality.append(f(population[-1]))
    best_indiv = copy.deepcopy(x0)
    try:
        for i in range(generations):
            for targetindex in range(len(population)):
                [v1,v2,v3] = random.sample(population,3)
                trial_vector = v1+v2-v3
                target_vector = population[targetindex]
                crossover(trial_vector,target_vector)
                bound(trial_vector)
                trial_quality = f(trial_vector)
                if trial_quality < quality[targetindex]:
                    population[targetindex] = trial_vector
                    quality[targetindex] = trial_quality
                
            # replace None results
            for j in range(len(population)):
                while quality[j] == None:
                    population[j] = new_individual()
                    quality[j] = f(population[j])

            best_index = quality.index(min(quality))
            if disp == 1:
                print "generation",str(i+1).rjust(7),"     f =",quality[best_index]
                print population[best_index]
            best_indiv = copy.deepcopy(population[best_index])
    except KeyboardInterrupt:
        if disp == 1:
            print "Stopping computation"
    return best_indiv



if __name__ == "__main__":
    #def f(x):
    #    return (x[0] - 4)**2 + (x[1]-8)**2 + (x[2]-1)**2 + (x[3]-9)**2 + (x[4]-8)**2 + (x[5]-2)**2
    
    #fmin_differential_evolution(f,scipy.array([1,1,1,1,1,1]),10000,20000,[[1e-1,1e1],[1e-1,1e1],[1e-1,1e1],[1e-1,1e1],[1e-1,1e1],[1e-1,1e1]],0.2)


    #x0 = numpy.array(medians)
    #x0 = scipy.array([0,0,0,0,0,0])
    #boundaries = [(0,10),(0,10),(0,10),(0,10),(0,10),(0,10)]
    #islogs     = [False,False,False,False,False,False]
    fmin_gen(f_opt,somedians,population_size=10,survivors=20,generations=2,bounds=boundaries,variable_is_logarithmic=islogs,intruders=0,use_pp=True)


