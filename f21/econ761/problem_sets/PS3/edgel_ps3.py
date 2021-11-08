# -*- coding: utf-8 -*-
"""
    This file performs all analayses requested in part 2 of problem set 3
    for Econ 761
    
    Date created:  01 Nov 2021
    Last modified: 02 Nov 2021
    Author: Danny Edgel
"""

# Load the pyBLP package and any other module necessary for conducting the
# analysis
import pyblp as blp 
import numpy as np
import pandas as pd

# save base directory for I/O
baseDir = 'C:/Users/edgel/Google Drive/UW-Madison/f21/econ761/problem_sets/PS3'

# load the data
ps3     = pd.read_excel((baseDir + '/data/cereal_ps3.xls'))
demog   = pd.read_excel((baseDir + '/data/demog_ps3.xls'))
    
v       = pd.read_csv((baseDir + '/data/v.csv'), header = None, 
                      index_col = False,
                      names = [('c' + str(i)) for i in range(1, 80)])

# transform demographic data to prepare it for estimation
demog = pd.concat([demog, v.reset_index(drop = True)], axis = 1)

demog['weights'] = 1/20
varnames = ['income', 'income_sq', 'age', 'child']
newnames = {}
vnames   = {}
for j in range(0, len(varnames)):
    varname = varnames[j]
    agt = 0
    for i in range(j*20 + 1, (j+1)*20 + 1):
        agt += 1
        newnames[('v' + str(i))] = (varname + str(agt))
        vnames[('c' + str(i))] = ('v' + str(j) + 'a' + str(agt))
        
newstubs = ['v0a', 'v1a', 'v2a', 'v3a']
for s in newstubs:
    varnames.append(s)        
demog = demog.rename(columns = newnames)
demog = demog.rename(columns = vnames)
demog = pd.wide_to_long(demog, i = ['city', 'year', 'quarter'], j = 'agent',
                        stubnames = varnames).reset_index()

nodes = {'v0a': 'nodes0', 'v1a': 'nodes1', 'v2a': 'nodes2', 'v3a': 'nodes3'}
demog.rename(columns = nodes, inplace = True)

# rename key variables and generate market_ids to make pyBLP happy
ps3.rename(columns = {'price': 'prices', 'share': 'shares',
                      'firm': 'firm_ids'}, inplace = True)

ps3['market_ids'] = ('C' + 
                     (ps3['city'].astype(str)).apply(
                         lambda x: '{0:0>2}'.format(x)) + 'Q' + 
                     ps3['quarter'].astype(str))

demog['market_ids'] = ('C' + 
                       (demog['city'].astype(str)).apply(
                           lambda x: '{0:0>2}'.format(x)) + 'Q' + 
                       demog['quarter'].astype(str))

znames = {}
for i in range(1, 20):
    znames[('z' + str(i))] = ('demand_instruments' + str(i - 1))
    
ps3.rename(columns = znames, inplace = True)

# define formulations
x1_form = blp.Formulation('0 + prices', absorb='C(brand)')
x2_form = blp.Formulation('1 + prices + sugar + mushy')
forms   = (x1_form, x2_form)

agent_form = blp.Formulation('0 + income + income_sq + age + child')

# simulate 20 individuals

"""
    There appear to be some irregularities in the 'v' and 'demog' data from
    Nevo--they don't appear to represent what the document says they do;
    load the Nevo data directly to get sensible results
"""
demog = pd.read_csv(blp.data.NEVO_AGENTS_LOCATION)
agent_form = blp.Formulation('0 + income + income_squared + age + child')

problem = blp.Problem(forms, ps3, agent_form, demog)


# make initial guess for sigma and pi
sigma0 = np.diag([0.3302, 2.4526, 0.0163, 0.2441])
pi0 = np.array([
  [ 5.4819,  0,      0.2037,  0     ],
  [15.8935, -1.2000, 0,       2.6342],
  [-0.2506,  0,      0.0511,  0     ],
  [ 1.2650,  0,     -0.8091,  0     ]
])


# specify optimization options
bfgs = blp.Optimization('bfgs', {'gtol': 1e-5})

# solve the model
blp_results = problem.solve(
    sigma0,
    pi0,
    optimization = bfgs,
    method='1s'
)
blp_results


# calculate markups and marginal costs
mc      = blp_results.compute_costs()
markups = blp_results.compute_markups(costs = mc)





