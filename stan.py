from datetime import datetime
from copy import copy,deepcopy
import sys

import numpy as np
from scipy.stats import bernoulli
import matplotlib.pyplot as plt
import pystan
from sklearn.grid_search import GridSearchCV
from sklearn.cross_validation import ShuffleSplit
from sklearn.base import BaseEstimator
from sklearn.metrics import roc_curve,auc

from trillion import load_components,load_odorants_tests_results

with open('trillion.stan') as f:
    code = f.read()

class StanModel2(BaseEstimator):
    def __init__(self,**kwargs):
        for key,value in kwargs.items():
            setattr(self,key,value)

    def set_model(self,code):
        self.model = pystan.StanModel(model_code=code)
        #def __deepcopy__(self):
        #    return self
        #self.model.__deepcopy__ = __deepcopy__

    def set_data(self,data):
        self.data = data
        for key,value in data.items():
            setattr(self,key,value)

    '''
    def get_model(self):
        if type(self.model) is dict:
            print self.model
            print self.model.values()[1][0]
            print self.model.values()[1][0].values()[1][0]       
            print self.model.values()[1][0].values()[1][0].values()[1][0]          
            for value in self.model.values():
                if type(value) is list:
                    if value[0].__class__.__name__ == 'StanModel':
                        return value[0]
        else:
            return self.model
    '''

    def optimize(self,X,y):
        for key in self.data.keys():
            self.data[key] = getattr(self,key)
        self.data.update({'n_obs':X.shape[0],
                          'subject_ids':X[:,0],
                          'test_ids':X[:,1],
                          'correct':y})
        self.best = self.model.optimizing(data=self.data)

    def get_params(self,deep=False):
        return self.__dict__
    
    def fit(self,X,y):
        self.optimize(X,y)
        for key,value in self.best.items():
            setattr(self,key,value)
        
    def transform(self,X,y=None,**fit_params):
        return X

    def predict(self,X):
        """X should be an Nx2 array, with subject ids in the first column
        and test ids in the second column.
        """  
        subjects = X[:,0] 
        tests = X[:,1]
        n_samples = X.shape[0]
        prediction = np.zeros(n_samples)
        for i in range(n_samples):
            difficulty = self.r_ta[tests[i]-1]
            discriminability = self.r_a[tests[i]-1]
            skill = self.r_s[subjects[i]-1]
            prediction[i] = 0.333 + 0.667/(1+np.exp(-discriminability*(skill-difficulty)));
        return prediction

    def score(self,X,y):
        prediction = self.predict(X)
        return np.log(bernoulli.pmf(y,prediction)).sum()
    
    def hist(self):
        p_correct_hit = self.p_correct[np.where(self.correct==1)]
        p_correct_miss = self.p_correct[np.where(self.correct==0)]
        plt.hist(p_correct_hit)
        plt.hist(p_correct_miss)
        plt.show()

    def roc(self):
        fpr, tpr, _ = roc_curve(self.correct,self.p_correct)
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr,tpr)
        plt.show()
        return roc_auc

THEN = None
def tic():
    global THEN
    THEN = datetime.now()

def toc(activity='Something'):
    now = datetime.now()
    delta = now - THEN
    seconds = delta.days*24*3600 + delta.seconds + delta.microseconds/1e6
    print('%s took %.3g seconds' % (activity,seconds))

def get_posterior_means(fit):
    means = {}
    x = fit.extract()
    for key,value in x.items()[:-1]:
        means[key] = value.mean(axis=0)
    return means

def make_data(search_data=None,test_size=0):
    components = load_components()
    odorants,tests,results = load_odorants_tests_results(components)
    n_subjects = len(set([x.subject_id for x in results]))
    n_tests = len(tests)
    n_obs = int(len(results)*(1-test_size))
    n_molecules = 128
    
    mixtures1 = np.zeros((n_tests,n_molecules))
    mixtures2 = np.zeros((n_tests,n_molecules))
    for test_id in tests:
        single = [components.index(molecule) for molecule in tests[test_id].single.components]
        for molecule in single:
            mixtures1[test_id-1,molecule] = 1
        double = [components.index(molecule) for molecule in tests[test_id].double.components]
        for molecule in double:
            mixtures2[test_id-1,molecule] = 1
    tests = np.array([result.test.id for result in results])
    subjects = np.array([result.subject_id for result in results])
    correct = np.array([int(result.correct) for result in results])
    
    data = locals()
    for to_del in ['components','odorants','results','single','double',
                   'test_id','molecule']:
        del data[to_del]
    if 'search_data' in data:
        del data['search_data']
    del data['test_size']
    if search_data:
        data.update({key:value[0] for key,value in search_data.items()})
    return data

if __name__ == '__main__':
    estimator = StanModel2()
    estimator.set_model(code)
    
    search_data = {'sigma_s':[0.3,1.0,3.0]}
    test_size = 0.1
    data = make_data(search_data=search_data,test_size=test_size)
    estimator.set_data(data)
    cv = ShuffleSplit(data['n_obs'],test_size=test_size)
    grid = GridSearchCV(estimator,search_data,cv=cv)
    
    y = data['correct']
    X = np.vstack((data['subject_ids'],data['test_ids'])).transpose()
    grid.fit(X,y)
        
    '''
    tic()
    model = pystan.StanModel(model_code=code)
    toc('Model compilation')
    data = make_data()    
    best = model.optimizing(data=data)
    fit = model.sampling(data=data)
    means = fit
    '''
