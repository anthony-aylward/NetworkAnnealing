#---------------------------------------------------------------------------------------#
#                                network_annealing.py                                   #
#---------------------------------------------------------------------------------------#

# Parameter optimization for network_sequence_simulator by simulated annealing. 

# This script assumes that the directory ~/HIVClustering/ is a clone of the eponymous git 
# repository. It outputs its results in the file ~/NetworkAnnealing/annealing_resultst.txt.

# First import a few dependencies

import argparse
import subprocess
import os.path
import csv
import numpy




#---------------------------------------------------------------------------------------#
#                                     Arguments                                         #
#---------------------------------------------------------------------------------------#

# Initialize the argument parser

arguments = argparse.ArgumentParser( 
            description='optimize transmission network parameters by simulated annealing')




#----------------------- arguments for input and output files ---------------------------#

arguments.add_argument('-s','--sequences', required=True,
        help='Provide the MSA with sequences which were used to make the distance file.')
arguments.add_argument('-f','--fasta', required=True,
        help='Destination path for output FASTA.')
arguments.add_argument('-t','--tn93', required=True,
        help='Destination path for output .csv.')




#-------------------------- arguments for network parameters ----------------------------#

# During the annealing process, each parameter will either be free to move or be held 
# fixed. This choice is specified by how the parameters are entered as arguments. For
# example, --size 100 indicates a network of fixed size 100, while --freesize 100 
# indicates an initial network size of 100 that will change during the annealing process.

# size

sizeMEG = arguments.add_mutually_exclusive_group()
sizeMEG.add_argument('-n','--size', type=int,
          help='Set value for a FIXED network size.')
sizeMEG.add_argument('--freesize', type=int,
          help='Set initial value for a FREE network size.')
          
# mean number of days between transmissions

daysMEG = arguments.add_mutually_exclusive_group()
daysMEG.add_argument('-d','--days', type=float,
          help='Set value for a FIXED mean transmission interval (days).')
daysMEG.add_argument('--freedays', type=float,
          help='Set initial value for a FREE mean transmission interval (days).')
          
# nucleotide substitution rate

rateMEG = arguments.add_mutually_exclusive_group()
rateMEG.add_argument('-r','--rate', type=float,
          help='Set value for a FIXED nucleotide substitution rate.')
rateMEG.add_argument('--freerate', type=float,
          help='Set initial value for a FREE nucleotide substitution rate.')
              
# number of starting lineages

linMEG = arguments.add_mutually_exclusive_group()
linMEG.add_argument('-l','--lineages', type=int,
          help='Set value for a FIXED number of starting lineages.')
linMEG.add_argument('--freelineages', type=int,
          help='Set initial value for a FREE number of starting lineages.')
          
# lineage initiation rate

splitMEG = arguments.add_mutually_exclusive_group()
splitMEG.add_argument('-p','--split', type=float,
          help='Set value for a FIXED lineage initiation rate.')
splitMEG.add_argument('--freesplit', type=float,
          help='Set initial value for a FREE lineage initiation rate.')

# random attachment rate

randMEG = arguments.add_mutually_exclusive_group()
randMEG.add_argument('-m','--random', type=float,
          help='Set value for a FIXED random attachment rate.')
randMEG.add_argument('--freerandom', type=float,
          help='Set initial value for a FREE random attachment rate.')
          
# subset size

subMEG = arguments.add_mutually_exclusive_group()
subMEG.add_argument('-u','--subset', type=int,
          help='Set value for a FIXED subset size.')
subMEG.add_argument('--freesubset', type=int,
          help='Set initial value for a FREE subset size.')
              
# attachment bias

biasMEG = arguments.add_mutually_exclusive_group()
biasMEG.add_argument('-b','--bias', type=float,
          help='Set value for FIXED attachment bias.')
biasMEG.add_argument('--freebias', type=float,
          help='Set initial value for FREE attachment bias.')
    
# sampling delay

sampMEG = arguments.add_mutually_exclusive_group()
sampMEG.add_argument('-y','--sampling', type=int,
          help='Set value for FIXED sampling delay (days).')
sampMEG.add_argument('--freesampling', type=int,
          help='Set initial value for FREE sampling delay (days).')
              
# burst size

burstMEG = arguments.add_mutually_exclusive_group()
burstMEG.add_argument('-x','--burst', type=float,
          help='Set value for FIXED mean burst size.')
burstMEG.add_argument('--freeburst', type=float,
          help='Set initial value for FREE mean burst size.')



#------------------------------------- parsing -----------------------------------------#

settings = arguments.parse_args()




#-------------------------- initializing param and freevars -----------------------------#

# Here we generate the dictionary param, which stores the inputs obtained from the 
# argument parser in string form, and the dictionary freevars, which stores the free
# variables in integer or floating point form.

param = {}
freevars = {}

param['sequences'] = settings.sequences
param['fasta'] = settings.fasta
param['tn93'] = settings.tn93

if settings.freesize or settings.freesize==0:
    freevars['size'] = int(settings.freesize)
    param['size'] = str(settings.freesize)
else:
    param['size'] = str(settings.size)
    
if settings.freedays or settings.freedays==0:
    freevars['days'] = float(settings.freedays)
    param['days'] = str(settings.freedays)
else:
    param['days'] = str(settings.days)
    
if settings.freerate or settings.freerate==0:
    freevars['rate'] = float(settings.freerate)
    param['rate'] = str(settings.freerate)
else:
    param['rate'] = str(settings.rate)
    
if settings.freelineages or settings.freelineages==0:
    freevars['lineages'] = int(settings.freelineages)
    param['lineages'] = str(settings.freelineages)
else:
    param['lineages'] = str(settings.lineages)
    
if settings.freesplit or settings.freesplit==0:
    freevars['split'] = float(settings.freesplit)
    param['split'] = str(settings.freesplit)
else:
    param['split'] = str(settings.split)
    
if settings.freerandom or settings.freerandom==0:
    freevars['random'] = float(settings.freerandom)
    param['random'] = str(settings.freerandom)
else:
    param['random'] = str(settings.random)
    
if settings.freesubset or settings.freesubset==0:
    freevars['subset'] = int(settings.freesubset)
    param['subset'] = str(settings.freesubset)
else:
    param['subset'] = str(settings.subset)
    
if settings.freebias or settings.freebias==0:
    freevars['bias'] = float(settings.freebias)
    param['bias'] = str(settings.freebias)
else:
    param['bias'] = str(settings.bias)
    
if settings.freesampling or settings.freesampling==0:
    freevars['sampling'] = int(settings.freesampling)
    param['sampling'] = str(settings.freesampling)
else:
    param['sampling'] = str(settings.sampling)
    
if settings.freeburst or settings.freeburst==0:
    freevars['burst'] = float(settings.freeburst)
    param['burst'] = str(settings.freeburst)
else:
    param['burst'] = str(settings.burst)




#---------------------------------------------------------------------------------------#
#                               Simulating the network                                  #
#---------------------------------------------------------------------------------------#

# In this section, we call the network generator network_sequence_simulator using the 
# subprocess module and the parameters identified in the settings.

def simulate_network(param):         
    subprocess.call([
                     '/opt/python-3.3.1/bin/python3',
                     os.path.expanduser('~/HIVClustering/bin/network_sequence_simulator'),
                     '-s',param['sequences'],
                     '-f',param['fasta'],
                     '-t',param['tn93'],
                     '-n',param['size'],
                     '-d',param['days'],
                     '-r',param['rate'],
                     '-l',param['lineages'],
                     '-p',param['split'],
                     '-m',param['random'],
                     '-u',param['subset'],
                     '-b',param['bias'],
                     '-y',param['sampling'],
                     '-x',param['burst']
                    ])
                          
    # In case tn93 fails to run during network_sequence_simulator's run, we call it here.
                          
    subprocess.check_call([
                           '/usr/local/bin/tn93', 
                           '-q', 
                           '-t', str(0.015), 
                           '-o', param['tn93'], 
                           param['fasta']
                          ])
    
    # We'll also need to run tn93 again without threshold in order to compute the mean
    # distance in the overall network later.   
                          
    subprocess.check_call([
                           '/usr/local/bin/tn93', 
                           '-q', 
                           '-t', str(1), 
                           '-o', param['tn93'][0:-4]+'_nothresh.csv', 
                           param['fasta']
                          ])
    return



#---------------------------------------------------------------------------------------#
#                             Computing network metrics                                 #
#---------------------------------------------------------------------------------------#

# Once the output .csv file has been generated by tn93, we can read it and use it to 
# compute the summarizing metrics for the inferred network.

def compute_metrics(filename):
    csv1 = open(filename, 'r')
    csv2 = open(filename[0:-4]+'_nothresh.csv', 'r')
    network = [edge for edge in csv.reader(csv1)][1:]
    networknothresh = [edge for edge in csv.reader(csv2)][1:]

    tails = [edge[0] for edge in network]
    heads = [edge[1] for edge in network]
    distances = [float(edge[2]) for edge in networknothresh]
    nodes = list(set(tails+heads))
    degrees = [(tails+heads).count(node) for node in nodes]

    metrics = {}
    metrics['nodes'] = len(nodes)
    metrics['edges'] = len(tails)
    metrics['meandist'] = sum(distances)/len(distances)
    if sum([numpy.log(degree) for degree in degrees]) == 0:
        metrics['charexp'] = float('inf')
    else:
        metrics['charexp'] = 1 + len(nodes)/sum([numpy.log(degree) for degree in degrees])
    
    csv1.close()
    csv2.close()
    
    return metrics




#---------------------------------------------------------------------------------------#
#                            Computing objective function                               #
#---------------------------------------------------------------------------------------#

# With the network's metrics available, we can compare them to the target metrics and 
# compute the value of our objective function (sum of squared errors on each normalized 
# metric).

def compute_objective(metrics):
    objnodes = 34
    objedges = 54
    objdist = 0.0583
    objchar = 3
    objective = sum([
                     ((metrics['nodes']-objnodes)/objnodes)**2,
                     ((metrics['edges']-objedges)/objedges)**2,
                     ((metrics['meandist']-objdist)/objdist)**2,
                     ((metrics['charexp']-objchar)/objchar)**2
                    ])
    return objective



#---------------------------------------------------------------------------------------#
#                                       Annealing                                       #
#---------------------------------------------------------------------------------------#

# Here we set the cooling schedule for the annealing process and construct the annealing
# step on which we will iterate later. First we construct a linear cooling schedule.

coolingrate = 0.1
schedule = [coolingrate*x for x in range(1,751)][::-1]

# At each annealing step, each free parameter will be updated with a new value drawn from
# a distribution with mode equal to the parameter's current value and a fixed variance
# proportional to the parameter's initial value. The variance for each free parameter is
# defined below.

def set_variance(freevars):
    variancefactor = 0.1
    variance = {}
    for key in freevars.keys():
        variance[key] = freevars[key] * variancefactor
    return variance

# Now we define how each free variable should be updated during the annealing step. We 
# compute distributional parameters from the desired mode and variance and draw from a
# Gamma, Beta, binomial, negative binomial, or poisson distribution with those parameters.

def update_free_variables(freevars, variance):
    for key in freevars.keys():
        k = 1 + (freevars[key]**2)/(2*variance[key]) + numpy.sqrt((freevars[key]**2)/variance[key]+(freevars[key]**4)/(4*variance[key]**2))
        theta = 1 / (freevars[key]/(2*variance[key]) + numpy.sqrt(1/variance[key]+(freevars[key]**2)/(4*variance[key]**2)))
        if key in ['size','lineages','subset','sampling']:
            if theta < 1:
                p = 1-theta
                n = int(numpy.floor(k*theta/p))
                freevars[key] = int(max(1,numpy.random.binomial(n,p)))
            elif theta > 1:
                p = 1-1/theta
                n = int(numpy.ceil(k/p))
                freevars[key] = int(max(1,numpy.random.negative_binomial(n,p)))
            elif theta == 1:
                lam = k*theta
                freevars[key] = int(max(1,numpy.random.poisson(lam)))
        elif key in ['rate','split']:
            mode = 10*freevars[key]
            var = 10*variance[key]
            
            a = var * mode**(-3)
            b = var * (7-3*mode**(-1))*mode**(-2) - mode**(-1) + 1
            c = var * (7-2*mode**(-1))*(2-mode**(-1))*mode**(-1) + mode**(-1) - 2
            d = var * (2-mode**(-1))**2*(3-mode**(-1))
            
            alpha = max([1.1]+[float(x) for x in list(numpy.roots((a,b,c,d))) if numpy.imag(x)==0])
            beta = max([1.1, (1/mode - 1) * alpha - 1/mode + 2])
            
            freevars[key] = 1/10*numpy.random.beta(alpha,beta)
        else:
            freevars[key] = numpy.random.gamma(k,theta)
    return freevars

# Once the new free variables have been determined, param should be updated appropriately
# for use in the next call of network_sequence_simulator

def update_parameters(param):
    for key in param.keys():
        if key in freevars.keys():
            param[key] = str(freevars[key])
    return param

# Now we combine the above functions into a fully formed annealing step. We store the 
# current state, update the free variables, update param, and generate a new network
# using the new parameters. Then we compute metrics, compute the objective function, and
# make a decision about whether to keep the new state or return to the old state. Once 
# the decision has been made, the annealing step is complete and is ready to iterate 
# again.

def annealing_step(param,freevars,variance,objective,T):
    oldparam = dict(param)
    oldfreevars = dict(freevars)
    oldobj = float(objective)
    freevars = update_free_variables(freevars, variance)
    param = update_parameters(param)
    simulate_network(param)
    metrics = compute_metrics(param['tn93'])
    objective = compute_objective(metrics)
    probaccept = min(1, numpy.exp((oldobj - objective)/T))
    randomfloat = numpy.random.random()
    
    if randomfloat > probaccept:
        param = dict(oldparam)
        freevars = dict(oldfreevars)
        objective = float(oldobj)
    
    return param, freevars, objective




#---------------------------------------------------------------------------------------#
#                                       Iterating                                       #
#---------------------------------------------------------------------------------------#

# Here we define the overall iterating structure of the script. We set the variances for
# the free variables, perform an initial annealing step to initialize the objective value,
# and then begin iteration along the temperature schedule.

variance = set_variance(freevars)
objective = float('inf')

T = schedule[0]
bestparam, bestfreevars, bestobjective = annealing_step(param,freevars,variance ,objective,T)

for temp in schedule:
    T = float(temp)
    param, freevars, objective = annealing_step(param,freevars,variance,objective,T)
    print(freevars)
    if objective < bestobjective:
        bestparam, bestfreevars, bestobjective = param, freevars, objective
        
        
        

#---------------------------------------------------------------------------------------#
#                                        Output                                         #
#---------------------------------------------------------------------------------------#

# Once a 'best' set of parameters and objective have been determined, we generate and 
# write a .txt file that contains them.

outputfile = open(os.path.expanduser('~/NetworkAnnealing/annealing_results.txt'), 'w')

keys = ['size','days','rate','lineages','split','random','subset','bias','sampling','burst']

output = ''.join([
                  'ANNEALING RESULTS\n\n',
                  'PARAMETERS\n',
                  '\n'.join([key + ': ' + bestparam[key] for key in keys]),
                  '\n\nOBJECTIVE\n' + str(bestobjective) + '\n'
                 ])
       
outputfile.write(output)
outputfile.close()
