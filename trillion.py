from pylab import *
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
        if(self.N and components != self.N):
            pass #print "Odorant does not contain %d components." % self.N
        self.components = components if components else []
        self.vector = array([i for i in range(self.C) if i in self.components])

    C = 10 # Number of candidate molecules that can be selected from
          # to build the odorant.  
    N = 3 # Number of molecules that an odorant must contain.  

    def r(self,other):
        if len(self.components) == len(other.components):
            return self.hamming(other)/2
        else:
            return None

    def hamming(self, other):
        diff = sum(abs(self.vector - other.vector))
        return diff

    def add_component(self, component):
        self.components.append(component)
        self.vector[component] = 1

    def remove_component(self, component):
        self.components.remove(component)
        self.vector[component] = 0

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
    def __init__(self,test_uid,odorants,dilution,correct):
        self.id = test_uid
        self.odorants = odorants
        self.dilution = dilution
        self.correct = correct

    def add_odorant(self,odorant):
        self.odorants.append(odorant)

class Result(object):
    def __init__(self,test,subject_id,correct):
        self.test = test
        self.subject_id = subject_id
        self.correct = correct

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
    for row in reader:
        uid,n,r,percent,dilution,correct = row[:6]
        component_names = [x for x in row[6:36] if len(x)]
        component_names = [x.replace('4-Methyl-3-penten-2-one',
                           '4-methylpent-3-en-2-one') for x in component_names]
        outcomes = row[36:62]  
        if uid.isdigit():
            dilution = DILUTION[dilution]
            #correct = CORRECT[correct]
            odorant_key = hash(tuple(component_names))
            if odorant_key not in odorants:
                components = [component for component in all_components \
                              if component.name in component_names]
                if len(components) not in [1,10,20,30]:
                    print uid,[x for x in component_names if x not in [y.name for y in components]]
                odorant = Odorant(components)
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
    return odorants,tests,results #.values(),tests.values()

def main():
    components = load_components()
    odorants,tests = load_odorants_tests_results(components)
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

if __name__ == '__main__':
    history = main()
    maxx,meen,empty = zip(*history)
    #plot(minn)
