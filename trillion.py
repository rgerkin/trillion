

import sys
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.misc import factorial
from scipy.stats import binom
import random
import csv
from datetime import datetime
from collections import OrderedDict,Counter

DILUTION = {'1/4':0.25, '1/2':0.5, 'not diluted':1.0}
CORRECT = {'right':True, 'wrong':False}
ALPHAS = [0.0005,0.001,0.002,0.005,0.01,0.02,0.05]
N_REPLICATES_LIST = np.arange(2,16)**2
C_LIST = 2**np.arange(5,15)
VERBOSE = False

def print_(*args,**kwargs):
    global VERBOSE
    if VERBOSE:
        print(*args,**kwargs)

def list_combos(n,k):
    """
    Returns a list of all combinations of n choose k.
    """
    
    return list(combinations(list(range(n)),k))

def get_n_combos(n,k):
    """
    Returns the number of combinations of n choose k.
    Uses Stirling's approximation when numbers are very large.
    """

    result = None
    try:
        fac_n = factorial(n)
        fac_k = factorial(k)
        fac_nk = factorial(n-k)
        summ = fac_n + fac_k + fac_nk
        if summ == np.inf or summ != summ:
            raise ValueError("Values too large. Using Stirling's approximation")
        result = fac_n/fac_k
        result /= fac_nk
    except: # Catch all large number exceptions.  
        pass
    if not result or result == np.inf or result != result: # No result yet.  
        x = stirling(n) - stirling(k) - stirling(n-k) # Use Stirling's approx.
        result = np.exp(x) 
    
    return result

def stirling(n):
    """
    Given n, returns Stirling's approximation for log(n!).
    """

    return n*np.log(n) - n + 0.5*np.log(2*np.pi*n)

def sphere(N,C,R):
    """
    Formula for sphere from Bushdid supplemental material.
    N = Number of components in a mixture.  
    C = Number of components to choose from.  
    R = Number of differing components.  
    """

    if R == 0:
        result = 1
    else:
        result = get_n_combos(N,min(R,N))*get_n_combos(C-N,min(R,C-N))
    return result

def ball(N,C,R):
    """
    Formula for ball from Bushdid supplemental material.
    N = Number of components in a mixture.  
    C = Number of components to choose from.  
    R = Maximum number of differing components.  
    """

    result = 0
    for r in range(0,R+1):
        result += sphere(N,C,r)
    return result

def disc(N,C,d):
    """
    Formula for number of discriminable odors from Bushdid supplemental material.
    N = Number of components in a mixture.  
    C = Number of components to choose from.  
    d = Discriminability limen.  
    """

    low = get_n_combos(C,N) / ball(N,C,int(np.floor(d/2)))
    high = get_n_combos(C,N) / ball(N,C,int(np.ceil(d/2)))
    result = low+(high-low)*(d/2-np.floor(d/2)) # Interpolate.  
    if result < 1:
        result = 1
    return result

def fdr(alpha,p_list):
    """
    Controls for false discovery rate using the Benjamin and Hochberg procedure.
    Given a nominal Type 1 error rate alpha and a list of nominal p-values,
    returns a corresponding list of booleans, set to True only if the null 
    hypothesis should still be rejected after controlling the false discovery rate. 
    """
    
    m = len(p_list)
    p_list_sorted = sorted(p_list)
    reject = [False for p in p_list]
    for k,p in enumerate(p_list):
        print_(p,k*alpha/m)
        if p <= k*alpha/m:
            reject[p_list.index(p)] = True
    return reject

class Odorant(object):
    """
    A mixture of molecules, defined by the presence of absence of the 
    candidate molecules in the mixture.
    """

    def __init__(self, components=None):
        """
        Builds odorant from a list of components.
        """
        
        self.components = components if components else []
        
    name = None # Name of odorant, built from a hash of component names.  

    C = 128 # Number of components from which to choose.  

    @property
    def N(self):
        """
        Number of components in this odorant.
        """

        return len(self.components)

    def r(self,other):
        """
        Number of replacements (swaps) to get from self to another other odorant.
        """

        if len(self.components) == len(other.components):
            return self.hamming(other)/2
        else:
            return None

    def overlap(self,other,percent=False):
        """
        Overlap between self and another odorant.  Complement of r. 
        Optionally report result as percent relative to number of components.
        """

        overlap = self.N - self.r(other)
        if percent:
            overlap = overlap*100.0/self.N
        return overlap

    def hamming(self, other):
        """
        Hamming distance between self and another odorant.  
        Synonymous with number of d, the number of total 'moves' to go from
        one odorant to another.
        """

        x = set(self.components)
        y = set(other.components)
        diff = len(x)+len(y)-2*len(x.intersection(y))
        return diff

    def add_component(self, component):
        """
        Adds one component to an odorant.
        """
        
        self.components.append(component)
        
    def remove_component(self, component):
        """
        Removes one component to an odorant.
        """
        
        self.components.remove(component)
        
    def descriptor_list(self,source):
        """
        Given a data source, returns a list of descriptors about this odorant.
        """
        
        descriptors = []
        for component in self.components:
            if source in component.descriptors:
                desc = component.descriptors[source]
                if type(desc) == list:
                    descriptors += desc
                if type(desc) == dict:
                    descriptors += [key for key,value in list(desc.items()) if value > 0.0]
        return list(set(descriptors)) # Remove duplicates.  

    def descriptor_vector(self,source,all_descriptors):
        """
        Given a data source, returns a vector of descriptors about this odorant.
        The vector will contain positive floats.
        """
        
        vector = np.zeros(len(all_descriptors))
        for component in self.components:
            if source in component.descriptors:
                desc = component.descriptors[source]
                if type(desc) == list:
                    for descriptor in desc:
                        index = all_descriptors.index(descriptor)
                        assert index >= 0
                        vector[index] += 1
                if type(desc) == dict:
                    this_vector = np.array([value for key,value in sorted(desc.items())])
                    vector += this_vector
        return vector

    def descriptor_vector2(self,all_descriptors):
        """
        Returns a vector of descriptors about this odorant, combining multiple
        data sources.  
        """
        
        n_descriptors_dravnieks = len(all_descriptors['dravnieks'])
        n_descriptors_sigma_ff = len(all_descriptors['sigma_ff'])
        vector = np.zeros(n_descriptors_dravnieks+n_descriptors_sigma_ff)
        for component in self.components:
            if 'dravnieks' in component.descriptors:
                desc = component.descriptors['dravnieks']
                this_vector = np.array([value for key,value in sorted(desc.items())])
                vector[0:n_descriptors_dravnieks] += this_vector
            elif 'sigma_ff' in component.descriptors:
                desc = component.descriptors['sigma_ff']
                for descriptor in desc:
                    index = all_descriptors['sigma_ff'].index(descriptor)
                    assert index >= 0
                    vector[n_descriptors_dravnieks+index] += 1
        return vector
    
    def described_components(self,source):
        """
        Given a data source, returns a list of the components which are
        described by that source, i.e. those that have descriptors.
        """
        
        return [component for component in self.components \
                if source in component.descriptors]

    def n_described_components(self,source):
        """
        Given a data source, returns the number of components that are
        described by that data source.
        """

        return len(self.described_components(source))

    def fraction_components_described(self,source):
        """
        Given a data source, returns the fraction of components that are
        described by that data source.
        """

        return self.n_described_components(source) / self.N

    def __str__(self):
        """
        String representation of the odorant.
        """
        
        return ','.join([str(x) for x in self.components])

class Component(object):
    """
    A single molecule, which may or may not be present in an odorant.
    """

    def __init__(self,component_id,name,cas,percent,solvent):
        """
        Components are defined by a component_id from the Bushdid et al
        supplemental material, a name, a CAS number, a percent dilution, 
        and a solvent.
        """

        self.id = component_id
        self.name = name
        self.cas = cas
        self.percent = percent
        self.solvent = solvent
        self.descriptors = {} # An empty dictionary.  

    def set_descriptors(self,source,cas_descriptors):
        """
        Given a data source, sets descriptors for this odorant using
        a dictionary where CAS numbers are keys, and descriptors are values.
        """

        assert type(source)==str and len(source)
        if self.cas in cas_descriptors:
            self.descriptors[source] = cas_descriptors[self.cas]
            # For sigma_ff this will be a list.  
            # For dravnieks this will be a dict.  

class Test(object):
    """
    One kind of experimental test performed by Bushdid et al.  
    This corresponding to a 'triangle test' with two odorants, and is
    defined by those odorants.
    """
    
    def __init__(self,test_uid,odorants,dilution,correct):
        """
        Tests are defined by their universal identifier (UID), the 3 odorants 
        used (2 should be identical), the dilution, and the identity of the
        correct response, which should be the odd-ball.  
        """

        self.id = test_uid
        self.odorants = odorants
        self.dilution = dilution
        self.correct = correct

    def add_odorant(self,odorant):
        """
        Adds one odorant to this test.
        """
        
        self.odorants.append(odorant)

    @property
    def double(self):
        """
        Returns the odorant present twice in this test.
        """
        
        for odorant in self.odorants:
            if self.odorants.count(odorant) == 2:
                return odorant
        return None

    @property
    def single(self):
        """
        Returns the odorant present once in this test.
        """
        
        for odorant in self.odorants:
            if self.odorants.count(odorant) == 1:
                return odorant
        return None

    @property
    def pair(self):
        """
        Returns the odorant pair in this test, with the odorant present
        twice listed first.
        """
        
        return (self.double,self.single)

    @property
    def N(self):
        """
        Returns the number of components in each of the odorants.  
        This a single value since they should all have the same number
        of components.
        """

        return self.double.N

    @property
    def r(self):
        """
        Returns the number of component replacements (swaps) separating one of
        the odorants from the other.
        """

        return self.double.r(self.single)

    def overlap(self, percent=False):
        """
        Returns the overlap (complement of r) between the two odorants.
        Optionally returns this as a percentage of N.
        """

        return self.double.overlap(self.single,percent=percent)

    @property
    def common_components(self):
        """
        Returns a list of components common to the two odorants.
        """

        d = set(self.double.components)
        s = set(self.single.components)
        return list(s.intersection(d))

    @property
    def unique_components(self):
        """
        Returns a list of components that exactly one of the two odorants has.
        """

        d = set(self.double.components)
        s = set(self.single.components)
        return list(s.symmetric_difference(d))

    def unique_descriptors(self,source):
        """
        Given a data source, returns a list of descriptors that 
        exactly one of the two odorants has.
        """

        sl = self.single.descriptor_list(source)
        dl = self.double.descriptor_list(source)
        unique = set(dl).symmetric_difference(set(sl))
        return list(unique)

    def common_descriptors(self,source):
        """
        Given a data source, returns a list of descriptors that 
        are common to the two odorants.
        """

        sl = self.single.descriptor_list(source)
        dl = self.double.descriptor_list(source)
        unique = set(dl).intersection(set(sl))
        return list(unique)

    def descriptors_correlation(self,source,all_descriptors):
        """
        Given a data source, returns the correlation between the descriptors
        of the two odorants.
        """

        sv = self.single.descriptor_vector(source,all_descriptors)
        dv = self.double.descriptor_vector(source,all_descriptors)
        return np.corrcoef(sv,dv)[1][0]

    def descriptors_correlation2(self,all_descriptors):
        """
        Returns the correlation between the descriptors
        of the two odorants, combining multiple data sources.
        """
        
        sv = self.single.descriptor_vector2(all_descriptors)
        dv = self.double.descriptor_vector2(all_descriptors)
        return np.corrcoef(sv,dv)[1][0]

    def n_undescribed(self,source):
        """
        Given a data source, returns the number of components from among the 
        two odorants that are not described by that source.
        """

        d = self.double.n_described_components(source)
        s = self.single.n_described_components(source)
        return (self.N-d,self.N-s)

class Result(object):
    """
    A test result, corresponding to one test given to one subject.
    """

    def __init__(self, test, subject_id, correct):
        """
        Results are defined by the test to which they correspond, 
        the id of the subject taking that test, and whether the subject 
        gave the correct answer.
        """ 
        
        self.test = test
        self.subject_id = subject_id
        self.correct = correct

class Distance(object):
    """
    An odorant distance, corresponding to distance between two odorants.
    No particular implementation for computing distance is mandated.
    """

    def __init__(self, odorant_i, odorant_j, distance):
        self.odorant_i = odorant_i
        self.odorant_j = odorant_j
        self.distance = distance

def load_components():
    """
    Loads all odorant components from Supplemental Table 1 of Bushdid et al.
    """

    components = []
    f = open('Bushdid-tableS1.csv','r',encoding='latin1')
    reader = csv.reader(f)
    next(reader)
    component_id = 0
    for row in reader:
        name,cas,percent,solvent = row[:4]
        if len(name):
            component = Component(component_id,name,cas,percent,solvent)
            components.append(component)
            component_id += 1
        else:
            break
    return components

def load_odorants_tests_results(all_components):
    """
    Given all odor components, loads the odorants, tests, and test results
    from Supplemental Table 2 of Bushdid et al.
    """
    odorants = {}
    tests = {}
    results = []
    f = open('Bushdid-tableS2.csv','r',encoding='latin1')
    reader = csv.reader(f)
    next(reader)
    row_num = 0
    for row in reader:
        uid,n,r,percent,dilution,correct = row[:6]
        component_names = [x for x in row[6:36] if len(x)]
        # The next line is required to account for inconsistent naming of one
        # of the components across the two supplemental tables.
        component_names = [x.replace('4-Methyl-3-penten-2-one',
                           '4-methylpent-3-en-2-one') for x in component_names]
        outcomes = row[36:62]  
        if uid.isdigit():
            uid = int(uid)
            dilution = DILUTION[dilution]
            odorant_key = hash(tuple(component_names))
            if odorant_key not in odorants:
                components = [component for component in all_components \
                              if component.name in component_names]
                if len(components) not in [1,10,20,30]:
                    # If an odorant has a number of components which is not
                    # either 1, 10, 20, or 30.  
                    print_(uid,[x for x in component_names if x not in [y.name for y in components]])
                odorant = Odorant(components)
                odorant.name = odorant_key
            elif row_num % 3 == 0:
                # If any component is repeated across all the tests.  
                print_("Repeat of this odorant: %d" % odorant_key)
            odorants[odorant_key] = odorant    
            if uid not in tests:
                tests[uid] = Test(uid,[],dilution,correct)
            test = tests[uid]
            test.add_odorant(odorant)
            if correct == 'right':
                test.correct = tests[uid].odorants.index(odorant)
        if len(outcomes[0]):
            for i,outcome in enumerate(outcomes):
                result = Result(test,i+1,CORRECT[outcome])
                results.append(result)
        row_num += 1
    return odorants,tests,results

def odorant_distances(results,subject_id=None):
    """
    Given the test results, returns a dictionary whose keys are odorant pairs
    and whose values are psychometric distances between those pairs, 
    defined as the fraction of discriminations that were incorrect.
    This can be limited to one subject indicated by subject_id, or else 
    by default it pools across all subjects.  
    """

    distances = {}
    distance_n_subjects = {} 
    for result in results:
        if subject_id and result.subject_id != subject_id:
            continue
        pair = result.test.pair
        if pair not in distances:
            distances[pair] = 0
            distance_n_subjects[pair] = 0
        distances[pair] += 0.0 if result.correct else 1.0
        distance_n_subjects[pair] += 1
    for pair in list(distances.keys()):
        # Divided by the total number of subjects.
        distances[pair] /= distance_n_subjects[pair]
    return distances

def ROC(results, N):
    """
    Given test results and a number of components N, returns a distribution 
    of the number of distinct components 'r' for correct trials (right) and 
    incorrect trials (wrong), in tests using odorants with N total components.
    These can later be plotted or used to generated an ROC curve.
    """
    
    right = []
    wrong = []
    for result in results:
        if result.test.N == N:
            r = result.test.r
            if result.correct:
                right.append(r)
            else:
                wrong.append(r)
    right = np.array(right) # Distribution of r for correct trials.  
    wrong = np.array(wrong) # Distribution of r for incorrect trials.  
    return (right,wrong)

def correct_matrix(results,N,overlap):
    """
    Given test results, a number of components N, and a level of overlap
    between odorants, returns a num_subjects by num_test matrix of booleans 
    corresponding to the correctness of that subject's response on that test.  
    """

    results = [r for r in results if r.test.N==N and r.test.overlap()==overlap]
    subjects = [r.subject_id for r in results]
    subjects = list(set(subjects))
    tests = [r.test for r in results]
    tests = list(set(tests))
    correct = np.zeros((len(subjects),len(tests)))
    correct -= 1 # Set to make sure every point gets set to 0 or 1 later.  
    for result in results:
        i = subjects.index(result.subject_id)
        j = tests.index(result.test)
        correct[i,j] = result.correct
    return correct

def fraction_disc(results,N,overlap,fig,alpha=None,multiple_correction=False,n_replicates=None):
    """
    Given test results, a number of components N, a level of overlap between
    odorants, a reference figure panel ('a' or 'b'), an optional choice of 
    significance threshold alpha, whether or not to do multiple comparisons
    correction (false discovery rate method), and an optional new number of 
    replicates (subjects or tests), returns an array containing either the 
    fraction of correct responses (if alpha is None) or whether or not that 
    fraction is significantly above chance (if alpha is a number).
    This function assists with generating variants of Figs. 2B, 2C, 3A, 
    and 3B in Bushdid et al.
    """

    assert fig in ['a','b']
    correct = correct_matrix(results,N,overlap)
    if fig == 'a':
        dim = 1
    elif fig == 'b':
        dim = 0
    fract_correct = np.mean(correct,dim)
    if alpha is not None:
        if not n_replicates:
            n_replicates = correct.shape[dim] # n_subjects or n_tests.
        ps = 1.0 - binom.cdf(fract_correct*n_replicates,n_replicates,1.0/3)
        if multiple_correction == 'bonferroni':
            alpha = alpha/len(ps)
        if multiple_correction == 'fdr':
            ps = np.array([p*len(ps)/(k+1) for k,p in enumerate(sorted(ps))])
        fract_sig = ps < alpha/2
        return fract_sig
    else:
        return fract_correct

def fig2x(results,fig='b',plot=True):
    """
    Given test results, a reference figure panel ('b' or 'c'), plots the data
    summary corresponding to Fig. 2 from Bushdid et al.
    """

    assert fig in ('b','c')
    overlap_dict = {10:[9,6,3,0],
               20:[19,15,10,5,0],
               30:[29,20,10,0]}
    
    f, axarr = plt.subplots(1, 3, sharey=True)
    for i,(N,overlaps) in enumerate(overlap_dict.items()):
        X = []
        Y = []
        for j,overlap in enumerate(overlaps):
            fract = fraction_disc(results,N,overlap,chr(ord(fig)-1),alpha=None)
            Y += list(fract)
            counts = Counter(fract)
            observed = {_:0 for _ in fract}
            for value in fract:
                inc = (observed[value] + 1)/(counts[value] + 1)-0.5
                X += [j+inc]
                observed[value] += 1

        axarr[i].scatter(X,Y)
        axarr[i].set_xticks([0,1,2,3,4])
        overlaps = np.array(overlaps)*100.0/N
        axarr[i].set_xticklabels(['%.2g'%_ for _ in overlaps])
    axarr[0].set_ylabel("% correct")
    axarr[1].set_xlabel("% mixture overlap")
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.4)
    #plt.show()

def fig3x(results,fig='a',alpha=0.05,multiple_correction=False,n_replicates=None,plot=True):
    """
    Given test results, a reference figure panel ('a' or 'b'), an optional 
    choice of significance threshold alpha, whether or not to do multiple comparisons
    correction (false discovery rate method), and an optional new number of 
    replicates (subjects or tests), returns the percent overlap of components 
    at which 50 percent of subjects/tests discriminate/can be discriminated 
    significantly above chance levels.  It does this by generating an analogue 
    of Fig. 3A or 3B from Bushdid et al, and computing the linear regression as 
    in those figures, and identifying the point at which the line intersects 
    with a horizontal line going through 50 percent on the ordinate.  
    By default, it plots the corresponding figure.  
    """

    tens_overlap = 100.0*np.array((9,6,3,0))/10.0
    twenties_overlap = 100.0*np.array([19,15,10,5,0])/20.0
    thirties_overlap = 100.0*np.array([29,20,10,0])/30.0
    
    def do(N,x,results):
        y = np.zeros(len(x))
        for i,overlap in enumerate(x):
            f = fraction_disc(results,N,int(overlap*N/100.0),fig,alpha=alpha,multiple_correction=multiple_correction,n_replicates=n_replicates)
            y[i] = np.mean(f)
        return y*100.0

    tens = do(10,tens_overlap,results)
    twenties = do(20,twenties_overlap,results)
    thirties = do(30,thirties_overlap,results)

    if plot:
        plt.scatter(tens_overlap,tens,s=20,c='b')
        plt.scatter(twenties_overlap,twenties,s=20,c='r')
        plt.scatter(thirties_overlap,thirties,s=20,c='g')
        plt.xlim(100,-1)
        plt.ylim(0,100)
    overlap = np.concatenate((tens_overlap,twenties_overlap,thirties_overlap))
    percent_disc = np.concatenate((tens,twenties,thirties))

    A = np.array([np.ones(len(overlap)), overlap])
    w = np.linalg.lstsq(A.T,percent_disc)[0] # obtaining the parameters

    xi = np.arange(0,100)
    line = w[0]+w[1]*xi # Regression line
    if plot:
        plt.plot(xi,line,'k-') # Plotting the line
        plt.xlabel('% mixture overlap')
        if fig == 'a':
            plt.ylabel('% subjects that can discriminate')
        elif fig == 'b':
            plt.ylabel('% mixtures that are discriminable')
    if w[0] == 0.0 and w[1] == 0.0:
        overlap = 0.0
    else:
        overlap = (50.0 - w[0])/w[1]
    print_(('50%% discrimination at %.3g%% overlap for alpha = '+('%g' if alpha else '%s')) \
        % (overlap,alpha))
    if overlap > 100.0:
        overlap = 100.0
    if overlap < 0.0:
        overlap = 0.0
    if plt:
        plt.plot([overlap,overlap],[0,100],'b--')
    return overlap

def overlap(results,fig='a',alphas=0.05*10.0**np.arange(-2,0.25,0.25),multiple_correction=False,n_replicates=None):
    """
    Given test results, a reference figure panel ('a' or 'b'), a range of 
    significance thresholds alpha, whether or not to do multiple comparisons
    correction (false discovery rate method), and an optional new number of 
    replicates (subjects or tests), plots the percent components overlap at 
    which 50 percent of subjects/tests discriminate/can be discriminated 
    significantly above chance levels, versus alpha.  
    """

    overlaps = []
    for alpha in alphas:
        overlap = fig3x(results,fig=fig,alpha=alpha,multiple_correction=multiple_correction,n_replicates=n_replicates,plot=False)
        overlaps.append(np.max([0,overlap]))
    if multiple_correction == 'bonferroni':
        color ='r'
    elif multiple_correction == 'fdr':
        color = 'g'
    elif not multiple_correction:
        color ='b'
    plt.scatter(alphas,overlaps,s=20,c=color)
    plt.xlim(0.1,np.min(alphas)*0.5)
    plt.ylim(-1,70)
    plt.xscale('log')
    plt.xlabel(r'Significance criterion $\alpha$')
    plt.ylabel('%% overlap for 50% discrimination')

def num_odors_vs_alpha(results,fig='a',alphas=0.05*10.0**np.arange(-2,0.25,0.25),multiple_correction=False):
    """
    Given test results, a reference figure panel ('a' or 'b'), a range of 
    significance thresholds alpha, and whether or not to do multiple comparisons
    correction (false discovery rate method), plots the number of odors implied
    by the equations in the supplemental material of Bushdid et al.  
    Uses the intermediate N = 20 case for convenience.  
    """

    N = 20
    n_odors_list = []
    for alpha in alphas:
        overlap = fig3x(results,fig=fig,alpha=alpha,multiple_correction=multiple_correction,plot=False)
        n_odors = disc(N,128,N*(100-overlap)/100)
        n_odors_list.append(float(n_odors))
    if multiple_correction == 'bonferroni':
        color ='g'
    elif multiple_correction == 'fdr':
        color = 'r'
    elif not multiple_correction:
        color ='b'
    
    plt.scatter(alphas,n_odors_list,s=20,c=color)
    plt.xlim(0.1,np.min(alphas)*0.5)
    plt.ylim(np.min(n_odors_list)*0.1,np.max(n_odors_list)*10)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'Significance criterion $\alpha$')
    plt.ylabel('Estimated number of odors')
    return (alphas,n_odors_list)

def num_odors_vs_replicates(results,fig='a',n_replicates_list=N_REPLICATES_LIST):
    """
    Given test results, a reference figure panel ('a' or 'b'), and an array of 
    new numbers of replicates (subjects or tests), plots the number of odors 
    implied by the equations in the supplemental material of Bushdid et al.  
    Uses N = 30 to illustrate the extreme values obtained.  
    Does not apply a correction for multiple comparisons for convenience, and 
    because that correction was not evident in Bushdid et al.  
    """

    N = 30
    n_odors_list = []
    for n_replicates in n_replicates_list:
        overlap = fig3x(results,fig=fig,alpha=0.05,multiple_correction=None,n_replicates=n_replicates,plot=False)
        overlap = np.min([overlap,100])
        n_odors = float(disc(N,128,N*(100-overlap)/100))
        print_(n_replicates,overlap,n_odors)
        n_odors_list.append(n_odors)
    
    if fig == 'a':
        color ='b'
    else:
        color ='g'
    
    plt.scatter(n_replicates_list,n_odors_list,s=20,c=color)
    plt.xlim(np.min(n_replicates_list)*0.5,np.max(n_replicates_list)*2)
    plt.ylim(np.min(n_odors_list)*0.1,np.max(n_odors_list)*10)
    plt.xscale('log')
    plt.yscale('log')
    if fig == 'a':
        plt.xlabel('Number of subjects')
    elif fig=='b':
        plt.xlabel('Number of tests')
    plt.ylabel('Estimated number of discriminable odors')
    return (n_replicates_list,n_odors_list)

def num_odors_vs_C(results,fig='a',Cs=C_LIST):
    """
    Given test results, a reference figure panel ('a' or 'b'), and an array of 
    component library sizes (the value C in Bushdid et al), plots the number of 
    odors implied by the equations in the supplemental material.  
    Uses N = 30 to illustrate the extreme values obtained.  
    Does not apply a correction for multiple comparisons for convenience, and 
    because that correction was not evident in Bushdid et al.  
    """

    N = 30
    n_odors_list = []
    overlap = fig3x(results,fig=fig,alpha=0.05,multiple_correction=False,plot=False)
    for C in Cs:
        n_odors = float(disc(N,C,N*(100-overlap)/100))
        print_(C,overlap,n_odors)
        n_odors_list.append(n_odors)
    
    if fig == 'a':
        color ='b'
    else:
        color ='g'
    
    plt.scatter(Cs,n_odors_list,s=20,c=color)
    plt.xlim(np.min(Cs)*0.5,np.max(Cs)*2)
    plt.ylim(np.min(n_odors_list)*0.1,np.max(n_odors_list)*10)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Number of components (C) \n in the molecular library")
    plt.ylabel(r'Estimated number of discriminable stimuli $\hat{O}$')
    return (Cs,n_odors_list)

def num_odors_vs_replicates_and_C(results,
               fig='a',
               alpha=0.05,
               n_replicates_list=N_REPLICATES_LIST,
               Cs=C_LIST):
    """
    Given test results, a reference figure panel ('a' or 'b'), an array of 
    new numbers of replicates (subjects or tests), and an array of 
    component library sizes (the value C in Bushdid et al), makes a heatmap of
    the number of odors implied by the equations in the supplemental material.  
    Uses the intermediate N = 20 case for convenience.  
    Applies a correction for multiple comparisons.  
    """

    n_replicates_list = np.rint(n_replicates_list).astype(int)

    N = 20
    n_odors_array = np.zeros((len(n_replicates_list),len(Cs)))
    for i,n_replicates in enumerate(n_replicates_list):
        overlap = fig3x(results,fig=fig,alpha=alpha,n_replicates=n_replicates,multiple_correction='fdr',plot=False)
        for j,C in enumerate(Cs):
            n_odors = float(disc(N,C,N*(100-overlap)/100))
            print_(C,n_replicates,overlap,n_odors)
            if np.isnan(n_odors):
                print_(N,C,overlap)
            n_odors_array[i,j] = n_odors
    
    return (n_replicates_list,Cs,n_odors_array)

def num_odors_vs_alpha_and_C(results,
               fig='a',
               alphas=ALPHAS,
               Cs=C_LIST):
    """
    Given test results, a reference figure panel ('a' or 'b'), an array of 
    alphas, and an array of component library sizes (the value C in Bushdid 
    et al), makes a heatmap of the number of odors implied by the equations 
    in the supplemental material. Uses the intermediate N = 20 case for 
    convenience.  Applies a correction for multiple comparisons.  
    """

    alphas = np.round(alphas,5)

    N = 20
    n_odors_array = np.zeros((len(alphas),len(Cs)))
    for i,alpha in enumerate(alphas):
        overlap = fig3x(results,fig=fig,alpha=alpha,multiple_correction='fdr',plot=False)
        for j,C in enumerate(Cs):
            n_odors = float(disc(N,C,N*(100-overlap)/100))
            print_(C,alpha,overlap,n_odors)
            if np.isnan(n_odors):
                print_(N,C,overlap)
            n_odors_array[i,j] = n_odors

    return (alphas,Cs,n_odors_array)

def num_odors_vs_alpha_and_replicates(results,
               fig='a',
               alphas=ALPHAS,
               n_replicates_list=N_REPLICATES_LIST):
    """
    Given test results, a reference figure panel ('a' or 'b'), an array of 
    alphas, and an array of new numbers of replicates (subjects or tests), 
    makes a heatmap of the number of odors implied by the equations 
    in the supplemental material. Uses the intermediate N = 20 case for 
    convenience.  Applies a correction for multiple comparisons.  
    """

    alphas = np.round(alphas,5)
    n_replicates_list = np.rint(n_replicates_list).astype(int)
    
    N = 20
    n_odors_array = np.zeros((len(alphas),len(n_replicates_list)))
    for i,alpha in enumerate(alphas):
        for j,n_replicates in enumerate(n_replicates_list):
            overlap = fig3x(results,fig=fig,alpha=alpha,n_replicates=n_replicates,multiple_correction='fdr',plot=False)
            n_odors = float(disc(N,128,N*(100-overlap)/100))
            print_(alpha,n_replicates,overlap,n_odors)
            if np.isnan(n_odors):
                print_(N,C,overlap)
            n_odors_array[i,j] = n_odors
    
    return (alphas,n_replicates_list,n_odors_array)

def get_results():
    components = load_components()
    odorants,tests,results = load_odorants_tests_results(components)
    return results

def write_data(data,name):
        """
        Given X,Y pairs (data) and a base name for the file, 
        writes that data to disk.
        """

        import csv
        
        with open('%s.csv' % name,'w') as f:
            transposed = []
            writer = csv.writer(f)
            for X,Y in data:
                transposed += [X]
                transposed += [Y]
            transposed = zip(*transposed)
            writer.writerows(transposed)

def figs():
    """
    Generates all figures and writes data points to disk.  
    """

    r = get_results()
            
    # Using the method of Fig. 3a, for varying values of alpha.    
    data = []
    data += [num_odors(r,fig='a')]
    data += [num_odors(r,fig='a',multiple_correction='fdr')]
    write_data(data,'a')
    plt.savefig('a_odors')
    plt.close()

    # Using the method of Fig. 3b, for varying values of alpha.    
    data = []
    data += [num_odors(r,fig='b')]
    data += [num_odors(r,fig='b',multiple_correction='fdr')]
    write_data(data,'b')
    plt.savefig('b_odors')
    plt.close()

    # Using the methods of Fig. 3a and 3b, for varying numbers of replicates.  
    data = []
    data += [num_odors2(r,fig='a')]
    data += [num_odors2(r,fig='b')]
    write_data(data,'2')
    plt.savefig('odors2')
    plt.close()

    # Using the methods of Fig. 3a and 3b, for varying values of C.    
    data = []
    data += [num_odors3(r,fig='a')]
    data += [num_odors3(r,fig='b')]
    write_data(data,'3')
    plt.savefig('odors3')
    plt.close()

def scratch():
    components = load_components()
    odorants,tests,results = load_odorants_tests_results(components)
    r = results
    kind = 'a'
    x1,y1,z1 = num_odors4(r,fig=kind)
    x2,y2,z2 = num_odors5(r,fig=kind)
    x3,y3,z3 = num_odors6(r,fig=kind)
    
    def heatmap(x,y,z,x_label,y_label):
        import matplotlib.pyplot as plt
        z = np.log10(z)
        fig, ax = plt.subplots()
        heatmap = ax.pcolor(z)

        # put the major ticks at the middle of each cell
        ax.set_xticks(np.arange(z.shape[1])+0.5, minor=False)
        ax.set_yticks(np.arange(z.shape[0])+0.5, minor=False)

        ax.set_xticklabels(y, minor=False)
        ax.set_yticklabels(x, minor=False)

        ax.set_xlabel(y_label)
        ax.set_ylabel(x_label)
        
        cbar = fig.colorbar(heatmap)
        ticks = [float(i.get_text()) for i in cbar.ax.get_yticklabels()]
        tick_labels = [str('$10^{%d}$' % tick) for tick in ticks]
        cbar.ax.set_yticklabels(tick_labels)
        cbar.ax.set_ylabel('# of discriminable stimuli')
     
    replicates_str = '# subjects' if kind=='a' else '# tests'
    alpha_str = 'Significance threshold $\\alpha$'
    c_str = 'Odorant panel size (C)'
    heatmap(x1,y1,z1,replicates_str,c_str) 
    heatmap(x2,y2,z2,alpha_str,c_str)   
    heatmap(x3,y3,z3,alpha_str,replicates_str)   
    plt.show()
    return [(x1,y1,z1),(x2,y2,z2),(x3,y3,z3)]