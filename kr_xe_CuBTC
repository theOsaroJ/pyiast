#data
pure_data= pd.read_csv('Cu-BTC Pure.csv')
P= pure_data['P (pa)'].dropna().to_numpy() #Pressure
kr_q= pure_data['q_kr (mol/kg)'].dropna().to_numpy() #Absolute Adsorption in Mol/Kg Framework
xe_q= pure_data['q_xe (mol/kg)'].dropna().to_numpy()
kr_m= pure_data['q_kr(molecules/unit cell)'].dropna().to_numpy() #Absolute Adsorption in molecules/Unit Cell
xe_m= pure_data['q_xe(molecules/unit cell)'].dropna().to_numpy()
kr_Q= pure_data['kr_heat (KJ/mol)'].dropna().to_numpy() #Heat of Adsorption
xe_Q= pure_data['xe_heat (KJ/mol)'].dropna().to_numpy()

#Plotting Adsorption in mol/kg Framework
plt.figure(1, figsize=(11,8))
plt.plot(P, kr_q, label='Krypton',color='r')
plt.scatter(P, kr_q,color='r')
plt.plot(P,xe_q,label='Xenon',color='b')
plt.scatter(P,xe_q,color='b')
plt.xlabel('Pressure(Pa)',weight='bold',fontsize=20)
plt.ylabel('Absolute Adsorption(mol/kg)',weight='bold',fontsize=20)
# plt.xscale('log')
# plt.yscale('log')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid()
plt.legend(loc='upper left')
plt.title('Absolute Adsorption of Pure Components Krypton and Xenon in Cu-BTC', weight='bold', fontsize=15)

plt.figure(2, figsize=(11,8))
plt.plot(P, kr_m, label='Krypton',color='r')
plt.scatter(P, kr_m,color='r')
plt.plot(P,xe_m,label='Xenon',color='b')
plt.scatter(P,xe_m,color='b')
# plt.xscale('log')
# plt.yscale('log')
plt.xlabel('Pressure(Pa)',weight='bold')
plt.ylabel('Absolute Adsorption(molecules/unit cell)',weight='bold')
plt.grid()
plt.legend(loc='upper left')
plt.title('Absolute Adsorption(molecules/unit cell) of Pure Components Krypton and Xenon in Cu-BTC', weight='bold', fontsize=12)

#IAST
#Isotherms

kr_interpolated= pyiast.InterpolatorIsotherm(pure_data, loading_key='q_kr (mol/kg)',pressure_key='P (pa)',
                                      fill_value=pure_data['q_kr (mol/kg)'].max())

xe_interpolated= pyiast.InterpolatorIsotherm(pure_data, loading_key='q_xe (mol/kg)', pressure_key='P (pa)',
                                             fill_value=pure_data['q_xe (mol/kg)'].max())

kr_langmuir=pyiast.ModelIsotherm(pure_data, loading_key='q_kr (mol/kg)',pressure_key='P (pa)',
                                 model='Langmuir',optimization_method='Powell')
# print(kr_langmuir.params)
xe_langmuir=pyiast.ModelIsotherm(pure_data, loading_key='q_xe (mol/kg)',pressure_key='P (pa)',
                                 model='Langmuir', param_guess={"M":25,"K":1.9e-11})

print(xe_langmuir.params)
#plotting
plt.figure(figsize=(11,8))
p_plot = np.logspace(np.log(pure_data['P (pa)'].min()),np.log(pure_data['P (pa)'].max()),num=100)
plt.scatter(P, kr_q, label='Krypton',marker='^',s=50,color='r',clip_on=False)
plt.plot(p_plot, kr_interpolated.loading(p_plot), color='y', linewidth=3,
         label='Linear interpolation', linestyle='--')
plt.plot(p_plot, kr_langmuir.loading(p_plot), color='k', linewidth=1,
         label='Langmuir Fits')

plt.scatter(P, xe_q, label='Xenon',marker='o',s=50,color='b',clip_on=False)
plt.plot(p_plot, xe_interpolated.loading(p_plot), color='y', linewidth=3
         , linestyle='--')
plt.plot(p_plot, xe_langmuir.loading(p_plot), color='k', linewidth=1)

plt.xlabel("Pressure (Pa)", weight='bold',fontsize=20)
plt.ylabel("Absolute Adsorption (mol/kg)", weight='bold',fontsize=20)
plt.xscale("log")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
# plt.yscale("log")
plt.legend(loc='upper left')
plt.tight_layout()
plt.title('Comparison of Adsorption Models for Krypton and Xenon', fontweight='bold', fontsize=12)
plt.grid()
plt.show()

mol_fraction= np.array([0.8])
total_pressure=np.array([1,10,100,1000,10000,100000,1000000,10000000,2e+7,3e+7,4e+7,5e+7,6e+7,7e+7,8e+7,9e+7,1e+8])
partial_pressure=[]

#Calculating the Partial Pressure based on Mol Fraction
for i in range(len(mol_fraction)):
    for j in range(len(total_pressure)):
        partial_pressure_eq= np.array([mol_fraction[i], 1-mol_fraction[i]])*total_pressure[j]
        partial_pressure.append(partial_pressure_eq)

# Using partial pressures to calculate Molal loading with pyIAST
q=[]
for i in range(len(partial_pressure)):
    pp=partial_pressure[i] 
    qq= pyiast.iast(pp, [kr_langmuir,xe_langmuir],warningoff=True)
    q.append(qq)

#krypton and xenon data
kr_q_iast= []
xe_q_iast=[]
for i in range(len(q)):
    kkr= q[i][0]
    xxe= q[i][1]
    kr_q_iast.append(kkr)
    xe_q_iast.append(xxe)

print('Krypton RMSE=',kr_langmuir.rmse)
print('Xenon RMSE=',xe_langmuir.rmse)
print(partial_pressure)

#comparing to GCMC MIXTURE DATA
#GCMC MIXTURE Data
mixture_data= pd.read_csv('Cu-BTC Mixture.csv')
mixture_data.head()
kr_q_mix= mixture_data['q_kr (mol/kg)'].dropna().to_numpy()
xe_q_mix= mixture_data['q_xe (mol/kg)'].dropna().to_numpy()
kr_q_mix_heat= mixture_data['kr_heat (KJ/mol)']
xe_q_mix_heat= mixture_data['xe_heat (KJ/mol)']
total_heat= mixture_data['Total_heat(KJ/mol)']
P= mixture_data['P (pa)'].dropna().to_numpy()


#Comparing IAST Values to that of the mixtures
plt.figure(1, figsize=(11,8))
plt.plot(P, kr_q_mix, label='Krypton GMC',color='r')
plt.scatter(P,kr_q_iast, label='Kr-IAST',color='r')
plt.plot(P, xe_q_mix, label='Xenon GCMC', color='b')
plt.scatter(P, xe_q_iast, label='Xe-IAST', color='b')
plt.legend()
plt.xlabel('Pressure  (Pa)',weight='bold', fontsize=20)
plt.ylabel('Absolute Loading(mol/kg)',weight='bold',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xscale('log')
plt.grid()
plt.title('Comparison of GCMC to IAST', weight='bold', fontsize=20)


#USING GP ON IAST
P= mixture_data['P (pa)'].dropna().to_numpy() 

#libraries
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RationalQuadratic,ConstantKernel as C, RBF, DotProduct, ExpSineSquared
from sklearn.metrics import r2_score, mean_squared_error
#Training Data set Randomly
rng= np.random.RandomState(1)
'''
Randomly picks values to be used as the training data set
'''
kr_q_iast= np.array(kr_q_iast) #Converting from list to array
xe_q_iast= np.array(xe_q_iast)
P= np.array(P)

n= len(kr_q_iast)
# for i in range(n):
training_indices= rng.choice(np.arange(kr_q_iast.size), size= 20, replace=False)
kr_train, xe_train = kr_q_iast[training_indices], xe_q_iast[training_indices]
P_train= P[training_indices]

#testing dataset

P_test= [value for value in P if value not in P_train]
P_test= np.array(P_test)

#Model for both Krypton and Xenon
kernel= C()*RationalQuadratic(length_scale=30, alpha=0.5,length_scale_bounds=(1e-15,1e15),alpha_bounds=(1e-15,1e15))
model= GaussianProcessRegressor(kernel=kernel,n_restarts_optimizer=4)
kr_train=kr_train.reshape(-1,1)
xe_train=xe_train.reshape(-1,1)
P_train= P_train.reshape(-1,1)
P_test= P_test.reshape(-1,1)
P= P.reshape(-1,1)

#Fitting model for Krypton
model.fit(P_train, kr_train)
kr_fit, σ_kr= model.predict(P, return_std=True)
print('RSquared for Krypton=',r2_score( kr_q_iast,kr_fit))
kr_mse= mean_squared_error(kr_q_iast,kr_fit)
kr_rmse= math.sqrt(kr_mse)
print('Rmse for Krypton=',kr_rmse)


# Fitting model for Xenon
model.fit(P_train, xe_train)
xe_fit, σ_xe= model.predict(P, return_std=True)
print('RSquared for Xenon=',r2_score( xe_q_iast,xe_fit))
xe_mse= mean_squared_error(xe_q_iast,xe_fit)
xe_rmse= math.sqrt(xe_mse)
print('Rmse for Xenon=',xe_rmse)


#Plotting
plt.figure(figsize=(11,8))
plt.plot(P, kr_q_iast, color='r', label='Krypton IAST')
plt.scatter(P, kr_fit,color='r', label='Krypton GP Fit')
plt.plot(P, xe_q_iast, color='b', label='Xenon IAST')
plt.scatter(P, xe_fit,color='b', label='Xenon GP Fit')
plt.legend(loc='best')
plt.xlabel('Pressure (Pa)', weight='bold', fontsize=20)
plt.ylabel('Absolute Adsorption(mol/kg)', weight='bold', fontsize=20)
plt.grid()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xscale('log')
# plt.title('IAST vs Gaussian')
plt.show()


print(len(P_train))
print(len(P_test))
print(len(P))
