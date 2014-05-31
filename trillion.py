from __future__ import division

import sys
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.misc import factorial
from scipy.stats import binom
import random
import csv
from datetime import datetime
from collections import OrderedDict

DILUTION = {'1/4':0.25, '1/2':0.5, 'not diluted':1.0}
CORRECT = {'right':True, 'wrong':False}

def list_combos(n,k):
    return list(combinations(range(n),k))

def get_n_combos(n,k):
    result = None
    try:
        result = (factorial(n)/(factorial(k))/factorial(n-k))
    except:
        pass
    if not result or result == np.inf or result != result:
        x = stirling(n) - stirling(k) - stirling(n-k)
        result = np.exp(x) 
    return result

def stirling(n):
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
        result = get_n_combos(N,R)*get_n_combos(C-N,R)
    return result

def ball(N,C,R):
    """Formula for ball from Bushdid supplemental material.
    N = Number of components in a mixture.  
    C = Number of components to choose from.  
    R = Maximum number of differing components.  
    """
    result = 0
    for r in range(0,R+1):
        result += sphere(N,C,r)
    return result

def disc(N,C,d):
    """Formula for number of discriminable odors from Bushdid supplemental material.
    N = Number of components in a mixture.  
    C = Number of components to choose from.  
    d = Discriminability limen.  
    """
    low = get_n_combos(C,N) / ball(N,C,int(np.floor(d/2)))
    high = get_n_combos(C,N) / ball(N,C,int(np.ceil(d/2)))
    result = low+(high-low)*(d/2-np.floor(d/2)) # Interpolate.  
    return result

def fdr(alpha,p_list):
    """Controls for false discovery rate using the Benjamin and Hochberg procedure."""
    m = len(p_list)
    p_list_sorted = sorted(p_list)
    reject = [False for p in p_list]
    for k,p in enumerate(p_list):
        print p,k*alpha/m
        if p <= k*alpha/m:
            reject[p_list.index(p)] = True
    return reject

'''
def ncube(n,k,d):
    combos = list(combinations(range(n),k))
    ref_combo = set(combos[0])
    for i,combo1 in enumerate(combos):
        intersection = ref_combo.intersection(set(combo1))
        diff = k - len(intersection)
        #print combo1,diff
    s = 0
    found = True
    while found:
        s += 1
        print "Trying s=%d" % s
        found = False
        meta_combos = list(combinations(combos,s))
        for meta_combo in meta_combos:
            diffs = [k-len(set(a).intersection(set(b))) for a in meta_combo for b in meta_combo if a is not b]
            #print diffs
            if all([x>=d for x in diffs]):
                found = True
                break
    print s-1
'''

'''
class Odor(object):
    """An odor, defined as that which cannot be discriminated from itself, 
    but can be discriminated from all other odors.  May contain 1 or 
    more odorants."""

    def __init__(self, odorants=None):
        self.odorants = odorants if odorants else []

    max_hamming = 2 # Maximum hamming distance of two odorants in an Odor.  

    def hammings(self):
        result = zeros((self.n_odorants,self.n_odorants))
        for i,odorant_i in enumerate(self.odorants):
            for j,odorant_j in enumerate(self.odorants):
                if i < j:
                    result[i][j] = odorant_i.hamming(odorant_j)
                elif i == j:
                    result[i][j] = 0
                else:
                    result[i][j] = result[j][i]
        self._hammings = result
        return result

    def hamming(self):
        h = self._hammings
        if len(h) > 1:
            maxx = h.max()
            meen = (sum(h)-sum(diag(h)))/(len(h)**2 - len(h))
            return (maxx,meen)
        else:
            return None

    @property
    def n_odorants(self):
        return len(self.odorants)

    def add_odorant(self, odorant):
        self.odorants.append(odorant)
    
    def remove_odorant(self, odorant):
        self.odorants.remove(odorant)

    def __str__(self):
        return u';'.join([str(x) for x in self.odorants])
'''
    
class Odorant(object):
    """A mixture of molecules, defined by the presence of absence of the 
    candidate molecules in the mixture."""

    def __init__(self, components=None):
        self.components = components if components else []
        #self.vector = np.array([i in self.components for i in range(self.C)])

    name = None # Name of odorant, built from a hash of component names.  

    C = 128 # Number of components from which to choose.  

    @property
    def N(self):
        """Number of components in this odorant."""
        return len(self.components)

    def r(self,other):
        if len(self.components) == len(other.components):
            return self.hamming(other)/2
        else:
            return None

    def overlap(self,other,percent=False):
        overlap = self.N - self.r(other)
        if percent:
            overlap = overlap*100.0/self.N
        return overlap

    def hamming(self, other):
        x = set(self.components)
        y = set(other.components)
        diff = len(x)+len(y)-2*len(x.intersection(y))
        return diff

    def add_component(self, component):
        self.components.append(component)
        #self.vector[component] = 1

    def remove_component(self, component):
        self.components.remove(component)
        #self.vector[component] = 0

    def descriptor_list(self,source):
        descriptors = []
        for component in self.components:
            if source in component.descriptors:
                desc = component.descriptors[source]
                if type(desc) == list:
                    descriptors += desc
                if type(desc) == dict:
                    descriptors += [key for key,value in desc.items() if value > 0.0]
        return list(set(descriptors)) # Remove duplicates.  

    def descriptor_vector(self,source,all_descriptors):
        #print "Getting %s descriptor vector" % source
        #tic = datetime.now()
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
        #print "%s" % (datetime.now()-tic)
        return vector

    def descriptor_vector2(self,all_descriptors):
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
        return [component for component in self.components \
                if source in component.descriptors]

    def n_described_components(self,source):
        return len(self.described_components(source))

    def fraction_components_described(self,source):
        return self.n_described_components(source) / self.N

    def __str__(self):
        return u','.join([str(x) for x in self.components])

class Component(object):
    def __init__(self,component_id,name,cas,percent,solvent):
        self.id = component_id
        self.name = name
        self.cas = cas
        self.percent = percent
        self.solvent = solvent
        self.descriptors = {}

    def set_descriptors(self,source,cas_descriptors):
        assert type(source)==str and len(source)
        if self.cas in cas_descriptors:
            self.descriptors[source] = cas_descriptors[self.cas]
            # For sigma_ff this will be a list.  
            # For dravnieks this will be a dict.  

class Test(object):
    """A test class, corresponding to a triangle test with two odorants."""
    def __init__(self,test_uid,odorants,dilution,correct):
        self.id = test_uid
        #assert len(odorants)==3 # There must be three stimuli.  
        self.odorants = odorants
        #assert not any([x is None for x in self.pair]) # Exactly two identical.
        self.dilution = dilution
        self.correct = correct

    def add_odorant(self,odorant):
        """Adds one odorant to this test."""
        self.odorants.append(odorant)

    @property
    def double(self):
        """Returns the odorant present twice in this test."""
        for odorant in self.odorants:
            if self.odorants.count(odorant) == 2:
                return odorant
        return None

    @property
    def single(self):
        """Returns the odorant present once in this test."""
        for odorant in self.odorants:
            if self.odorants.count(odorant) == 1:
                return odorant
        return None

    @property
    def pair(self):
        """Returns the odorant pair in this test, with the odorant present
        twice listed first."""
        return (self.double,self.single)

    @property
    def N(self):
        return self.double.N

    @property
    def r(self):
        return self.double.r(self.single)

    def overlap(self, percent=False):
        return self.double.overlap(self.single,percent=percent)

    @property
    def N(self):
        return self.double.N

    @property
    def r(self):
        return self.double.r(self.single)

    @property
    def common_components(self):
        d = set(self.double.components)
        s = set(self.single.components)
        return list(s.intersection(d))

    @property
    def unique_components(self):
        d = set(self.double.components)
        s = set(self.single.components)
        return list(s.symmetric_difference(d))

    def unique_descriptors(self,source):
        sl = self.single.descriptor_list(source)
        dl = self.double.descriptor_list(source)
        unique = set(dl).symmetric_difference(set(sl))
        return list(unique)

    def common_descriptors(self,source):
        sl = self.single.descriptor_list(source)
        dl = self.double.descriptor_list(source)
        unique = set(dl).intersection(set(sl))
        return list(unique)

    def descriptors_correlation(self,source,all_descriptors):
        sv = self.single.descriptor_vector(source,all_descriptors)
        dv = self.double.descriptor_vector(source,all_descriptors)
        return np.corrcoef(sv,dv)[1][0]

    def descriptors_correlation2(self,all_descriptors):
        sv = self.single.descriptor_vector2(all_descriptors)
        dv = self.double.descriptor_vector2(all_descriptors)
        return np.corrcoef(sv,dv)[1][0]

    def n_undescribed(self,source):
        d = self.double.n_described_components(source)
        s = self.single.n_described_components(source)
        return (self.N-d,self.N-s)

class Result(object):
    def __init__(self, test, subject_id, correct):
        self.test = test
        self.subject_id = subject_id
        self.correct = correct

class Distance(object):
    def __init__(self, odorant_i, odorant_j, distance):
        self.odorant_i = odorant_i
        self.odorant_j = odorant_j
        self.distance = distance

def load_components():
    components = []
    f = open('Bushdid-tableS1.csv','r')
    reader = csv.reader(f)
    reader.next()
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
    odorants = {}
    tests = {}
    results = []
    f = open('Bushdid-tableS2.csv','r')
    reader = csv.reader(f)
    reader.next()
    row_num = 0
    for row in reader:
        uid,n,r,percent,dilution,correct = row[:6]
        component_names = [x for x in row[6:36] if len(x)]
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
                    print uid,[x for x in component_names if x not in [y.name for y in components]]
                odorant = Odorant(components)
                odorant.name = odorant_key
            elif row_num % 3 == 0:
                print "Repeat of this odorant: %d" % odorant_key
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
    return odorants,tests,results #.values(),tests.values()

def odorant_distances(results,subject_id=None):
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
    for pair in distances.keys():
        distances[pair] /= distance_n_subjects[pair]
    return distances

def mds(odorants,results):
    distances = odorant_distances(results,subject_id=None)
    n_odorants = len(odorants)
    data = np.zeros((n_odorants,n_odorants),dtype='float64')
    mask = np.ones((n_odorants,n_odorants))
    for (double,single),distance in distances.items():
        i = odorants.keys().index(double.name)
        j = odorants.keys().index(single.name)
        data[i][j] = distance
        mask[i][j] = 0
        data[j][i] = distance
        mask[j][i] = 0
    distance_matrix = np.ma.MaskedArray(data,mask)

    from sklearn import manifold
    from sklearn.metrics import euclidean_distances
    seed = np.random.RandomState(seed=3)
    mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, 
                       random_state=seed, dissimilarity="precomputed", n_jobs=1)

    pos = mds.fit_transform(distance_matrix)
    
    '''
    nmds = manifold.MDS(n_components=2, metric=False, max_iter=3000, eps=1e-12,
                    dissimilarity="precomputed", random_state=seed, n_jobs=1,
                    n_init=1)
    
    npos = nmds.fit_transform(distance_matrix, init=pos)
    '''

    #fitted_distance_matrix = euclidean_distances(pos)
    #original = distance_matrix.reshape(len(pos)**2)
    #fitted = fitted_distance_matrix.reshape(len(pos)**2)
    #plot(pos[:,0],pos[:,1],'.')

    npos = pos
    #return original,fitted
    
    from matplotlib.collections import LineCollection
    plt.scatter(npos[:, 0], npos[:, 1], s=20, c='b')
    segments = []
    for (double,single),distance in distances.items():
        i = odorants.keys().index(double.name)
        j = odorants.keys().index(single.name)
        segments.append(((npos[i,:]),(npos[j,:])))
    print segments
    lc = LineCollection(segments)
    #lc.set_array(distance_matrix.flatten())
    #lc.set_linewidths(0.5 * np.ones(len(segments)))
    plt.gca().add_collection(lc)
    plt.show()

def ROC(results, N):
    """Return a distribution of number of distinct components for 
    correct trials (right) and incorrect trials (wrong)."""
    right = []
    wrong = []
    for result in results:
        if result.test.N == N:
            r = result.test.r
            if result.correct:
                right.append(r)
            else:
                wrong.append(r)
    right = np.array(right)
    wrong = np.array(wrong)
    return (right,wrong)

def correct_matrix(results,N,overlap):
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
    assert fig in ['a','b']
    correct = correct_matrix(results,N,overlap)
    if fig == 'a':
        dim = 1
    elif fig == 'b':
        dim = 0
    fract_correct = np.mean(correct,dim)
    if alpha is not None:
        if not n_replicates:
            n_replicates = correct.shape[dim] # n_subjects
        ps = 1.0 - binom.cdf(fract_correct*n_replicates,n_replicates,1.0/3)
        if multiple_correction == 'bonferroni':
            alpha = alpha/len(ps)
        if multiple_correction == 'fdr':
            #print 'Raw: '+str(sorted(ps))
            ps = np.array([p*len(ps)/(k+1) for k,p in enumerate(sorted(ps))])
            #print 'FDR: '+str(ps)
        fract_sig = ps < alpha/2
        return fract_sig
    else:
        return fract_correct

def fig3x(results,fig='a',alpha=0.05,multiple_correction=False,n_replicates=None,plot=True):
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

    # plotting the line
    xi = np.arange(0,100)
    line = w[0]+w[1]*xi # regression line
    if plot:
        plt.plot(xi,line,'k-')
        plt.xlabel('% mixture overlap')
        if fig == 'a':
            plt.ylabel('% subjects that can discriminate')
        elif fig == 'b':
            plt.ylabel('% mixtures that are discriminable')
    overlap = (50.0 - w[0])/w[1]
    print '50%% discrimination at %.3g%% overlap for alpha = %g' \
        % (overlap,alpha)
    return overlap

def overlap(results,fig='a',alphas=0.05*10.0**np.arange(-2,0.25,0.25),multiple_correction=False,n_replicates=None):
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
    plt.xlabel('Significance criterion alpha')
    plt.ylabel('%% overlap for 50% discrimination')

def num_odors(results,fig='a',alphas=0.05*10.0**np.arange(-2,0.25,0.25),multiple_correction=False):
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
    plt.xlabel('Significance criterion alpha')
    plt.ylabel('Estimated number of odors')
    return (alphas,n_odors_list)

def num_odors2(results,fig='a',n_replicates_list=2**np.arange(2,9,0.5)):
    N = 30
    n_odors_list = []
    for n_replicates in n_replicates_list:
        overlap = fig3x(results,fig=fig,alpha=0.05,multiple_correction=False,n_replicates=n_replicates,plot=False)
        overlap = np.min([overlap,100])
        n_odors = float(disc(N,128,N*(100-overlap)/100))
        print n_replicates,overlap,n_odors
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
    plt.xlabel('Number of replications (subjects or distinct mixtures)')
    plt.ylabel('Estimated number of discriminable odors')
    return (n_replicates_list,n_odors_list)

def num_odors3(results,fig='a',Cs=2**np.arange(6,15)):
    N = 30
    n_odors_list = []
    overlap = fig3x(results,fig=fig,alpha=0.05,multiple_correction=False,plot=False)
    for C in Cs:
        n_odors = float(disc(N,C,N*(100-overlap)/100))
        print C,overlap,n_odors
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
    plt.xlabel('Number of components (C) available to mix odorants')
    plt.ylabel('Estimated number of odors')
    return (Cs,n_odors_list)

def figs():
    components = load_components()
    odorants,tests,results = load_odorants_tests_results(components)
    r = results

    '''
    fig3x(r,fig='a')
    fig3x(r,fig='a',multiple_correction='fdr')
    plt.savefig('a_overlap')
    plt.close()

    fig3x(r,fig='b')
    fig3x(r,fig='b',multiple_correction='fdr')
    plt.savefig('b_overlap')
    plt.close()
    '''

    def write_data(data,name):
        import csv
        from itertools import izip
        with open('%s.csv' % name,'w') as f:
            transposed = []
            writer = csv.writer(f)
            for X,Y in data:
                transposed += [X]
                transposed += [Y]
            transposed = izip(*transposed)
            writer.writerows(transposed)
            
    data = []
    data += [num_odors(r,fig='a')]
    data += [num_odors(r,fig='a',multiple_correction='fdr')]
    write_data(data,'a')
    plt.savefig('a_odors')
    plt.close()

    data = []
    data += [num_odors(r,fig='b')]
    data += [num_odors(r,fig='b',multiple_correction='fdr')]
    write_data(data,'b')
    plt.savefig('b_odors')
    plt.close()

    data = []
    data += [num_odors2(r,fig='a')]
    data += [num_odors2(r,fig='b')]
    write_data(data,'2')
    plt.savefig('odors2')
    plt.close()

    data = []
    data += [num_odors3(r,fig='a')]
    data += [num_odors3(r,fig='b')]
    write_data(data,'3')
    plt.savefig('odors3')
    plt.close()
 
def fit(results, components):
    """Basic OLS model to predict test results."""
    X = np.zeros((len(results),26+len(components)*2))
    Y = np.zeros((len(results),1))
    for i,result in enumerate(results):
        Y[i] = result.correct
        x = np.zeros(len(components))
        test_components = result.test.double.components+result.test.single.components
        for component in test_components:
            try:
                index = components.index(component)
            except ValueError:
                print "Couldn't find %s" % component
                sys.exit(0)
            else:
                x[index] += 1
        for j in range(len(components)):
            X[i,26+2*j] = x[j]==1
            X[i,26+2*j+1] = x[j]==2
        X[i,result.subject_id] = 1             
    from scikits.statsmodels.api import OLS
    from sklearn import linear_model as lm
    clf = lm.Lasso(0.001)#alpha = alpha)
    results = clf.fit(X,Y)
    for i,beta in enumerate(results.coef_):
        result = round(1000*beta)
        if abs(result) > 0.1:
            print i,result
    print Y.shape,X.shape
    model = OLS(Y,X)
    results = model.fit()
    print results.summary()
    return Y,X

def fit2(results, tests, components):
    """Regularized model to predict averaged test results."""
    tests = [value for key,value in tests.items() if value.N > 1]
    X = np.zeros((len(tests),4+128+128))
    Y = np.zeros((len(tests),1))
    for test in tests:
        test.n_results = 0
    for result in results:
        if result.test not in tests:
            continue
        i = tests.index(result.test)
        test = tests[i]
        Y[i,0] += result.correct
        #Y[i,1] += not result.correct
        X[i,0] = 1
        X[i,1] = test.N == 20
        X[i,2] = test.N == 30
        X[i,3] = test.overlap(percent=True)/100.0
        for j,component in enumerate(components):
            if component in test.unique_components:
                X[i,j+4] += 1
            if component in test.common_components:
                X[i,j+4+128] += 1
        test.n_results += 1
    for i,test in enumerate(tests):
        Y[i] = (Y[i]+0.0)/test.n_results
    import statsmodels.api as sm
    #model = sm.GLM(Y,X,family=sm.families.Binomial())
    #results = model.fit()
    #print results.summary()
    '''
    from sklearn.linear_model import Lasso
    for alpha in 10.0 ** np.arange(-3,0,0.1):
        clf_l1_LR = Lasso(alpha=alpha)
        clf_l1_LR.fit(X, Y)
        coef_l1_LR = clf_l1_LR.coef_.ravel()
        coef_l1_LR[0:4] = 1
        X_sparse = np.delete(X,np.where(coef_l1_LR == 0.0)[0],1)
        if len(X_sparse):
            model = sm.OLS(Y,X_sparse)
            results = model.fit()
            r2 = 1-results.ssr/results.centered_tss
            adj_r2 = 1 - (results.nobs-1)/results.df_resid * (1-r2) 
            print 'Alpha = %f; Components = %d; R^2 = %f; Adj. R^2 = %f; AIC = %f; BIC = %f' % \
               (alpha,X_sparse.shape[1]-4,r2,adj_r2,
                results.aic,results.bic)
            #print results.summary()
    '''
    X_overlap_only = X[:,0:4]
    model = sm.OLS(Y,X_overlap_only)
    results = model.fit()
    print results.summary()
    X = []
    Y = []
    for i,test in enumerate(tests):
        for component in test.common_components:
            X.append(components.index(component))
        Y += [results.resid[i] for i in range(len(test.common_components))]
    '''
    Y_mean = [mean(array(Y)[where(array(X)==i)[0]]) for i in range(128)]
    indices = range(len(components))
    indices.sort(key=lambda x:Y_mean[x])
    X_sorted = [indices.index(i) for i in X]
    X_mean = [indices.index(i) for i in range(len(components))]
    scatter(X_sorted,Y,s=5)
    scatter(X_mean,Y_mean,s=50,c='r')
    '''
    return results,X,Y

def fit3(results, tests, components, all_descriptors):
    """Naive bayes prediction of test results."""
    results = [result for result in results if result.test.N > 1]
    n_subjects = 26
    n_descriptors_sigma_ff = len(all_descriptors['sigma_ff'])
    n_descriptors_dravnieks = len(all_descriptors['dravnieks'])
    n_descriptors = n_descriptors_sigma_ff + n_descriptors_dravnieks
    n_components = len(components)
    X = np.zeros((len(results),6+n_subjects))#+len(components)*2))
    Y = np.zeros((len(results),1))
    n_results = len(results)
    for i,result in enumerate(results):
        prog(i,n_results)
        test = result.test
        Y[i] = result.correct
        X[i,0] = test.N == 10
        X[i,1] = test.N == 20
        X[i,2] = test.N == 30
        #X[i,3] = test.overlap(percent=True)/100.0
        #X[i,4] = (test.overlap(percent=True)/100.0)**2
        #X[i,5] = (test.overlap(percent=True)/100.0)**3
        X[i,3] = test.descriptors_correlation2(all_descriptors)**2
        X[i,4] = test.descriptors_correlation2(all_descriptors)**4
        X[i,5] = test.descriptors_correlation2(all_descriptors)**6
        X[i,5+result.subject_id] = 1
        #for component in result.test.unique_components:
        #    X[i,3+n_subjects+components.index(component)]=1
        #for component in result.test.common_components:
        #    X[i,3+n_subjects+len(components)+components.index(component)]=1
    from sklearn.naive_bayes import BernoulliNB
    from sklearn.cross_validation import train_test_split
    from sklearn import metrics
    Y = np.squeeze(Y)
    clf = BernoulliNB()
    scores = []
    for i in range(25):
        X_train,X_test,Y_train,Y_test = train_test_split(X,Y,test_size=0.3)
        clf.fit(X_train,Y_train)
        Y_pred = clf.predict(X_test)
        scores.append(metrics.accuracy_score(Y_test,Y_pred))
    print "\nClassification Accuracy: %f" % np.mean(np.array(scores))
    return clf,X,Y

def fit4(results,tests,components,all_descriptors):
    tests = [value for key,value in tests.items() if value.N > 1]
    #sources = ['sigma_ff','dravnieks']
    n_descriptors_sigma_ff = len(all_descriptors['sigma_ff'])
    n_descriptors_dravnieks = len(all_descriptors['dravnieks'])
    n_descriptors = n_descriptors_sigma_ff + n_descriptors_dravnieks
    X = np.zeros((len(tests),9))#+len(sources)*1))
    Y = np.zeros((len(tests),1))
    for test in tests:
        test.n_results = 0
    n_results = len(results)
    for n_result,result in enumerate(results):
        prog(n_result,n_results)
        if result.test not in tests:
            continue
        i = tests.index(result.test)
        test = tests[i]
        #print test.id,test.N,test.n_undescribed('dravnieks')
        Y[i,0] += result.correct
        #Y[i,1] += not result.correct
        X[i,0] = 1
        X[i,1] = test.N == 20
        X[i,2] = test.N == 30
        X[i,3] = test.overlap(percent=True)/100.0
        X[i,4] = (test.overlap(percent=True)/100.0)**2
        X[i,5] = (test.overlap(percent=True)/100.0)**3
        X[i,6] = test.descriptors_correlation2(all_descriptors)**2
        X[i,7] = test.descriptors_correlation2(all_descriptors)**4
        X[i,8] = test.descriptors_correlation2(all_descriptors)**6
        #for j,source in enumerate(sources):
        #    offset = 2*j
            #X[i,4+j] = len(test.unique_descriptors(source))
            #X[i,5+j] = len(test.unique_descriptors(source))/test.N
            #X[i,4+j] = test.descriptors_correlation(source,all_descriptors[source])**4
            #if X[i,4+j] == 1.0:
            #    print result.test.id
            #X[i,7+j] = len(test.common_descriptors(source))
            #X[i,8+j] = len(test.common_descriptors(source))/test.N
        #for descriptor in test.unique_descriptors:
        #    index = all_descriptors.index(descriptor)
        #    X[i,9+index] = 1
        test.n_results += 1
    for i,test in enumerate(tests):
        Y[i] = (Y[i]+0.0)/test.n_results

    # Delete columns with only zeroes.
    X = np.delete(X,np.where(sum(X,0)==0)[0],axis=1)  

    import statsmodels.api as sm
    from sklearn import linear_model as lm
    from sklearn.cross_validation import train_test_split
    from sklearn import metrics
    Y = np.squeeze(Y)

    for col in range(1,X.shape[1]):
        #pass
        norm = np.linalg.norm(X[:,col])
        #print col,norm
        #X[:,col] /= norm

    model = sm.OLS(Y,X)
    results = model.fit()
    print results.summary()

    X_train,X_test,Y_train,Y_test = train_test_split(X,Y,test_size=0.3)
    #clf = BernoulliNB()
    #clf.fit(X_train,Y_train)
    
    for lasso in 10 ** np.arange(-6,-2,0.25):
        clf = lm.Lasso(lasso)#alpha = alpha)
        results = clf.fit(X_train,Y_train)
        print "Lasso parameter: %f" % lasso
        print "Out of sample R^2: %f" % np.corrcoef(Y_test,results.predict(X_test))[1][0]**2
    
    clf = lm.Lasso(1)#alpha = alpha)
    results = clf.fit(X_train,Y_train)
    for i,beta in enumerate(results.coef_):
        result = round(1000*beta)
        if abs(result) > 0.1:
            print i,result

    cols = np.where(np.array(results.coef_)==0.0)[0]
    if cols[0] == 0:
        cols = cols[1:]
    #X = np.delete(X,cols,axis=1)  
    
    X_train,X_test,Y_train,Y_test = train_test_split(X,Y,test_size=0.5)
    model = sm.OLS(Y_train,X_train)
    results = model.fit()
    print results.summary()

    #Y_test = Y_train
    #X_test = X_train
    #ss_res = sum((Y_test-results.predict(X_test))**2)
    #ss_tot = sum((Y_test-np.mean(Y_test))**2)
    #print ss_res,ss_tot
    #r2 = 1 - ss_res/ss_tot
    #print 
    print "Out of sample R^2: %f" % np.corrcoef(Y_test,results.predict(X_test))[1][0]**2
    
    return results,X,Y

def loess(x,y,frac=0.2,it=None,scatter=True):
    from statsmodels.nonparametric.smoothers_lowess import lowess 
    y = np.array(y)
    x = np.array(x)
    y = y[x.argsort()] # Sort y according to order of x.  
    x.sort() # Sort x in place.  
    if it is not None: # Helps if you are getting NaN's in the output.  
        d = lowess(y,x,frac=frac,it=it)
    else:
        d = lowess(y,x,frac=frac)
    return d
    #if scatter:
    #    plot(x,y,'.',x,d[:,1],'-')
    #else:
    #    plot(x,d[:,1],'-')

def tic():
    global T
    T = datetime.now()

def toc():
    global T
    print '%s' % (datetime.now() - T)

def prog(num,denom):
    fract = float(num)/denom
    hyphens = int(round(50*fract))
    spaces = int(round(50*(1-fract)))
    sys.stdout.write('\r%.2f%% [%s%s]' % (100*fract,'-'*hyphens,' '*spaces))
    sys.stdout.flush()        

def main():
    components = load_components()
    import sigma_ff,dravnieks
    sigma_ff_descriptors = sigma_ff.data
    dravnieks_descriptors = dravnieks.data
    for component in components:
        component.set_descriptors("sigma_ff",sigma_ff_descriptors)
        component.set_descriptors("dravnieks",dravnieks_descriptors)
    all_descriptors = {'sigma_ff':sigma_ff.descriptors,
                       'dravnieks':dravnieks.descriptors}
    odorants,tests,results = load_odorants_tests_results(components)
    return fit3(results,tests,components,all_descriptors)
    #result = mds(odorants,results)
    #return result

    """
    n_iterations = 1000
    n_odors = 5 # Initial number of odors.  
    n_odorants = 25
    combos = list_combos(Odorant.C,Odorant.N)
    component_lists = random.sample(combos,n_odorants)
    odorants = [Odorant(x) for x in component_lists]
    odors = [Odor() for x in range(n_odors)]
    for odorant in odorants:
        odor = random.choice(odors)
        odor.add_odorant(odorant)
    history = []
    for iteration in range(n_iterations):
        #print "===================="
        for odor in odors:
            #print "--------------------"
            #print odor
            if len(odor.odorants) <= 1:
                continue
            hammings = odor.hammings()
            sums = sum(hammings,0)
            if(len(sums)):
                most = odor.odorants[argmax(sums)] # Location of odor most dissimilar to others.  
                #print hammings,argmax(sums),most
                other_odor = random.choice(odors)
                #print most,odors.index(odor),odors.index(other_odor)
                odor.remove_odorant(most)
                other_odor.add_odorant(most)
        maxx,meen,empty = 0.,0.,0
        for odor in odors:
            try:
                x,y = odor.hamming()
                maxx += x
                meen += y
            except TypeError:
                empty += 1
        maxx /= (n_odors-empty)
        meen /= (n_odors-empty)
        history.append((maxx,meen,empty))
        print "After iteration %d: %f, %f, %d" % (iteration,maxx,meen,empty)
    return history
    """

if __name__ == '__main__':
    results,X,Y = main()
    #maxx,meen,empty = zip(*history)
    #plot(minn)
