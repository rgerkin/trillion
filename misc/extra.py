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