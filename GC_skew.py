from pysb import *
from pysb.simulator import BngSimulator
from pysb.util import alias_model_components
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.animation as animation
from scipy.stats import norm
from scipy.stats import poisson

plt.rcParams['animation.ffmpeg_path'] = '/opt/local/bin/ffmpeg'

n_genes = 100

# Poisson distribution of gene lengths

# Make a change

delta = 5
x = np.arange(20,81,delta)
mu = 50
y = poisson.pmf(x,mu)*delta*n_genes
y = np.array([int(round(yy)) for yy in y])

x_poi = []
y_poi = []
for i in range(len(x)):
    if y[i] != 0:
        x_poi.append(x[i])
        y_poi.append(y[i])

# Exponential distribution of gene lengths

delta = 10
x = np.arange(0,125+1,delta)
lam = 1/25
y = lam*np.exp(-lam*x)*delta
y /= np.sum(y)
y *= n_genes
y = np.array([int(round(yy)) for yy in y])
y[0] += n_genes - np.sum(y)
x += 25

x_exp = []
y_exp = []
for i in range(len(x)):
    if y[i] != 0:
        x_exp.append(x[i])
        y_exp.append(y[i])

# Log-normal distribution of gene lengths

delta = 500
x = np.arange(delta/2,4000,delta)
log10_x = np.log10(x)
mu = 2.85
sigma = 0.25
y = norm.pdf(log10_x, loc=mu, scale=sigma)
mean = 10**(mu+sigma**2/2)
sdev = np.sqrt(10**(2*mu+sigma**2)*(10**(sigma**2)-1))
y = y/sum(y)*n_genes
y = np.array([int(round(yy)) for yy in y])
x = np.array([int(round(xx)) for xx in x/mean*50])

x_log = []
y_log = []
for i in range(len(x)):
    if y[i] != 0:
        x_log.append(x[i])
        y_log.append(y[i])

# PLOT

# for (x,y) in [(x_poi,y_poi), (x_exp,y_exp), (x_log,y_log)]:
#        
#     fig, ax1 = plt.subplots()
#     print(x)
#     print(y)
#     print(sum(y))
#     dx = x[1]-x[0]
#     ax1.bar(x-0.2*dx, y, width=0.4*dx, color='b')
#     ax1.set_xlabel('Gene length')
#     ax1.set_ylabel('Number')
#     ax1.annotate('Total genes = %d' % int(sum(y)), (0.5,0.9), 
#                  xycoords='axes fraction', fontsize=14)
#     ax1.annotate('', xy=(0.25, 1.05), xytext=(0.45, 1.05),
#                 arrowprops=dict(facecolor='b', edgecolor='b'), 
#                 xycoords='axes fraction')
#          
#     ax2 = ax1.twinx()
#          
#     ax2.bar(x+0.2*dx, 10**np.linspace(0,1,len(x),endpoint=True), width=0.4*dx, color='g')
# #     ax2.set_yscale('log')
#     # ax2.xlabel('Gene length')
#     ax2.set_ylabel('k_skew')
#     ax2.annotate('', xy=(0.75, 1.05), xytext=(0.55, 1.05),
#                 arrowprops=dict(facecolor='g', edgecolor='g'), 
#                 xycoords='axes fraction')
#      
# plt.show()
# quit()

x = np.array(x_exp)
y = np.array(y_exp)

N_STATES = 2*x+1
# N_STATES = [101]*len(x)
N_GENES = y
log10_k_skew_min = 1
log10_k_skew_max = 1
K_SKEW = 10**np.linspace(log10_k_skew_min,log10_k_skew_max,len(x),endpoint=True)

def create_model(N_STATES,K_SKEW):
    Model()

    Monomer('G', ['dir', 'skew'], {'dir' : ['CD', 'HO'], 'skew' : ['_%d' % i for i in range(N_STATES)]})
    
    Parameter('Gtot', 1)
    Initial(G(dir='CD', skew='_%d' % ((N_STATES-1)/2)), Gtot)
    
    k_skew = [Parameter('kf_skew_HO_%d_%d' % (i,i+1), K_SKEW) for i in range(N_STATES-1)]
    [Rule('GC_skew_HO_%d_%d' % (i,i+1), G(dir='HO', skew='_%d' % i) >> G(dir='HO', skew='_%d' % (i+1)), k_skew[i]) for i in range(N_STATES-1)]
    
    Parameter('k_skew_CD', 10**log10_k_skew_min) 
    [Rule('GC_skew_CD_%d_%d' % (i,i+1), G(dir='CD', skew='_%d' % i) >> G(dir='CD', skew='_%d' % (i+1)), k_skew_CD) for i in range(N_STATES-1)]
    
    Parameter('k_cd_ho', 1)
    Parameter('k_ho_cd', 1.5)
    [Rule('Inversion_%d_%d' % (i,(N_STATES-1)-i), G(dir='CD', skew='_%d' % i) | G(dir='HO', skew='_%d' % ((N_STATES-1)-i)), k_cd_ho, k_ho_cd)  for i in range(N_STATES)]
    
    Observable('HO', G(dir='HO'))
    
    return model

t_end = 500
tspan = np.linspace(0, t_end, 5*t_end+1)
# sim = BngSimulator(model, tspan, verbose=True)

# n_levels = 10
# min_gcskew = 0.25
# x = np.arange(1,n_levels+1)
# 
# # linear: y = ax+b
# a = (min_gcskew-1)/(n_levels-1)
# b = 1-a
# y_lin = a*x+b
# 
# # exponential: y = exp(-kx)
# k = -np.log(min_gcskew)/(n_levels-1)
# y_exp = np.exp(-k*(x-1))
# 
# # y = 1/(1+c*(x-1))
# c = (1-min_gcskew)/min_gcskew/(n_levels-1)
# y_inv = 1/(1+c*(x-1))
# 
# y = y_exp
# k_remove = [int(round(z)) for z in (N_STATES-1)/2*(1-y)]
# 
# plt.bar(x, y)
# plt.yticks([0,0.25,0.5,0.75,1])
# plt.xlabel('batch')
# plt.ylabel('max GC skew')
 
# plt.show()
# quit()

n_runs = n_genes
len_tspan = len(tspan)
HO_tot = np.array(len_tspan*[0])
skew = np.empty((len_tspan, n_runs), dtype=float)
color = np.empty((len_tspan, n_runs), dtype=str)

# batch_size = int(n_runs/n_levels)
np.random.seed(44)
# for batch in range(n_levels):
start = 0
for batch in range(len(N_GENES)):
#     size = batch_size if batch < n_levels-1 else n_runs-(n_levels-1)*batch_size
    size = N_GENES[batch]
    model = create_model(N_STATES[batch],K_SKEW[batch])
    sim = BngSimulator(model, tspan, verbose=True)
    x = sim.run(n_runs=size, 
                netgen_args={'max_iter':1000000}, 
                seed=np.random.choice(int(1e6), size=size, replace=False) #,
#                 param_values={ k.name : 0 for k in kf_skew[len(kf_skew)-k_remove[batch]:] }
#                         if k_remove[batch] > 0 else None
                )
    observables = [x.observables] if size == 1 else x.observables
    species = [x.species] if size == 1 else x.species
    for t in range(len_tspan):
        for n in range(size):
            # Total HO genes
            HO_tot[t] += observables[n][t][0]
            # Skew and color 
            idx = np.where(species[n][t] == species[n][t].max())[0][0]
            sk = int(model.species[idx].monomer_patterns[0].site_conditions['skew'][1:])
#             skew[t][start+n] = (sk-((N_STATES[batch]-1)/2)) / ((N_STATES[batch]-1)/2)
            skew[t][start+n] = (sk-((N_STATES[batch]-1)/2)) / ((max(N_STATES)-1)/2)
            dir = model.species[idx].monomer_patterns[0].site_conditions['dir']
            color[t][start+n] = 'b' if dir == 'CD' else 'r'
    start += size

##### PLOTS #####

# Running average
# HO_tot = np.array([np.mean(HO_tot[i*int(len_tspan/10):(i+1)*int(len_tspan/10)]) for i in range(10)])
# t_mean = np.array([np.mean(tspan[i*int(len_tspan/10):(i+1)*int(len_tspan/10)]) for i in range(10)])
# plt.figure()
# plt.plot(t_mean, 1-HO_tot/n_runs, 'b', lw=2, label='Co-directional')
# plt.plot(t_mean, HO_tot/n_runs, 'r', lw=2, label='Head-on')
# plt.xlabel('Time average (a.u.)')
# plt.ylabel('Percent')
# plt.legend(loc=0)
# 
# # Average GC skew
# plt.figure()
# plt.plot(tspan, np.mean(skew, axis=1), 'k', lw=2)
# plt.xlabel('Time (a.u.)')
# plt.ylabel('Avg. GC skew')

# plt.show()
# quit()
#################

# Sort the GC skew values and colors
skew_sorted = np.empty_like(skew)
color_sorted = np.empty_like(color)
for t in range(len_tspan):
    idxs = np.flip(skew[t].argsort())
    skew_sorted[t] = skew[t][idxs]
    color_sorted[t] = color[t][idxs]

skew_sorted_CD = []
skew_sorted_HO = []
for t in range(len_tspan):
    skew_sorted_CD.append([])
    skew_sorted_HO.append([])
    for n in range(n_runs):
        if color_sorted[t][n] == 'b':
            skew_sorted_CD[t].append(skew_sorted[t][n])
        else:
            skew_sorted_HO[t].append(skew_sorted[t][n])
    skew_sorted_CD[t] += (n_runs-len(skew_sorted_CD[t]))*[0.]
    skew_sorted_HO[t] += (n_runs-len(skew_sorted_HO[t]))*[0.]

##### ANIMATION
fig = plt.figure() #figsize=(10,10))
ax = fig.add_subplot(111) #, aspect='equal')
ax.set_xlabel('gene')
ax.set_ylabel('GC skew')
# bars = ax.bar(range(n_runs), n_runs*[0], width=0.9, color='b')
bars1 = ax.bar(range(n_runs), n_runs*[0], width=0.8, color='b')
bars2 = ax.bar(range(n_runs), n_runs*[0], width=0.8, color='r')
ax.set_yticks(ticks=[-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1])
ax.set_ylim((-1.05,1.05))
cd = Rectangle((0, 0), 1, 1, fc="b")
ho = Rectangle((0, 0), 1, 1, fc="r")
plt.legend([cd,ho], ['Co-directional', 'Head-on'], loc=(0.7,1.05)) 
time_text = ax.text(0.1,1.05,'t = %g' % tspan[0], fontsize=25, color='r', transform=ax.transAxes)

def animate(t):
    print(tspan[t])
    time_text.set_text('t = %g' % tspan[t])
#     for i,bar in enumerate(bars):
#         bar.set_height(skew_sorted[t][i])
#         bar.set_color(color_sorted[t][i])
#         bar.set_height(skew[t][i])
#         bar.set_color(color[t][i])
    for i,b in enumerate(zip(bars1,bars2)):
        b[0].set_height(skew_sorted_CD[t][i])
        b[1].set_height(skew_sorted_HO[t][i])
    plt.tight_layout()
    return bars1+bars2

ani = animation.FuncAnimation(fig, animate, range(len(tspan)), interval=50, blit=True, repeat=False) 
ani.save('anim.mp4', dpi=100)

# plt.show()
