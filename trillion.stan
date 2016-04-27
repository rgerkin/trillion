functions {
        real p_disc(real dist, real kd, real exponent) { 
            return 0.333 + 0.666/(1+(kd/dist)^exponent);
            }
        }

data {
    int<lower=0,upper=100> n_dimensions; // number of perceptual dimensions
    int<lower=0,upper=200> n_molecules; // number of distinct molecules
    int<lower=0,upper=100> n_subjects; // number of subjects
    int<lower=0,upper=1000> n_tests; // number of tests
    int<lower=0> norm; // L1 or L2 norm for computing distance
    vector<lower=0,upper=1>[n_molecules] mixtures1[n_tests]; // Molecular presence or absence in the 1-sample mixtures 
    vector<lower=0,upper=1>[n_molecules] mixtures2[n_tests]; // Molecular presence or absence in the 2-sample mixtures 
    int<lower=1> n_obs;
    int<lower=1> subjects[n_obs]; // subject_id
    int<lower=1> tests[n_obs]; // test_id
    int<lower=0> correct[n_obs]; // 1 if correct, 0 if not
    }
    
parameters {
    vector<lower=0>[n_subjects] kd; // Subject kd's
    vector<lower=0>[n_subjects] exponent; // Subject exponents
    //vector<lower=0,upper=1>[n_subjects] p_dumb; // Subject dumbness probabilties
    matrix<lower=0,upper=1>[n_dimensions,n_molecules] molecule_coords;
    //int<lower=0,upper=1> molecule_bools[n_dimensions,n_molecules];
    
    real<lower=0> mu_kd;
    //real<lower=0> beta_dumb;
    real<lower=0> mu_exponent;
    //real<lower=0> molecule_alpha;
    //real<lower=0> molecule_beta;
    }
    
transformed parameters {
    vector<lower=0>[n_dimensions] mixture1_coords[n_tests];
    vector<lower=0>[n_dimensions] mixture2_coords[n_tests];
    real<lower=0> D[n_tests];
    real logp;
    
    for(i in 1:n_tests) {
        mixture1_coords[i] <- molecule_coords * mixtures1[i]; // compute mixture 1 coordinates
        mixture2_coords[i] <- molecule_coords * mixtures2[i]; // compute mixture 1 coordinates
        if(norm == 2) {
            D[i] <- distance(mixture1_coords[i],mixture2_coords[i]); // euclidean distance between mixtures
        }
        if(norm == 1) {
            D[i] <- 0;
            for(dim in 1:n_dimensions) {
                D[i] <- D[i] + fabs(mixture1_coords[i][dim] - mixture2_coords[i][dim]); // manhattan distance between mixtures
                }
            }
        }
    logp <- get_lp();    
    }

model {
    int s;
    int t;
    
    mu_kd ~ lognormal(0,100);
    mu_exponent ~ lognormal(0,100);
    //beta_dumb ~ uniform(1,100);
    //p_dim ~ beta(1,100);
    
    kd ~ lognormal(mu_kd,1.0);
    exponent ~ lognormal(mu_exponent,0.1);
    //p_dumb ~ beta(1,beta_dumb);

    for(dim in 1:n_dimensions) {
        molecule_coords[dim] ~ uniform(0,1);//exponential(molecule_alpha);//,molecule_beta);
        //molecule_coords[dim] ~ beta(molecule_alpha,molecule_beta);
    }
    /*
    for(mol in 1:n_molecules) {
        for(dim in 1:n_dimensions) {
            if(molecule_coords[dim,mol] == 0) {
                increment_log_prob(bernoulli_log(0,p_dim) +
                                   uniform_log(0,0,1));
                }
            else {
                increment_log_prob(bernoulli_log(1,p_dim) +
                                   uniform_log(molecule_coords[dim,mol],0,1));
                }
            }
        }
    */

    for(i in 1:n_obs) {
        s <- subjects[i];
        t <- tests[i];
        correct[i] ~ bernoulli(p_disc(D[t], kd[s], exponent[s]));
        }
    }

