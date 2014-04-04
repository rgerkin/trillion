import sys
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.misc import factorial
import random
import csv

DILUTION = {'1/4':0.25, '1/2':0.5, 'not diluted':1.0}
CORRECT = {'right':True, 'wrong':False}

def list_combos(n,k):
    return list(combinations(range(n),k))

def get_n_combos(n,k):
    return factorial(n)/(factorial(k)*factorial(n-k))

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

    def __str__(self):
        return u','.join([str(x) for x in self.components])

class Component(object):
    def __init__(self,component_id,name,cas,percent,solvent):
        self.id = component_id
        self.name = name
        self.cas = cas
        self.percent = percent
        self.solvent = solvent

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

class Result(object):
    def __init__(self,test,subject_id,correct):
        self.test = test
        self.subject_id = subject_id
        self.correct = correct

class Distance(object):
    def __init__(self,odorant_i,odorant_j,distance):
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

def ROC(results,N):
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

def fit(results,components):
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

def main():
    components = load_components()
    odorants,tests,results = load_odorants_tests_results(components)
    result = mds(odorants,results)
    return result

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
    history = main()
    #maxx,meen,empty = zip(*history)
    #plot(minn)
