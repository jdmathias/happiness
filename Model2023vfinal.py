########################################################################################
######################################## Happiness Treadmill model #####################
########################################################################################
import numpy as np
from scipy import fftpack
import matplotlib.pyplot as pl
pl.rcParams.update({'font.size': 22})
################################## Time parameters #####################################
days=10; hours = days*24; # Total time of the simulation
############################## Vector initialization ###################################
P=np.ones(hours)*0;# Total pleasure
A=np.ones(hours)*0;# Total afflictive affects
Ps=np.ones(hours)*0;# Pleasure provided (or not) by stimuli
As=np.ones(hours)*0;# Afflictive affects provided by stimuli
gammaP=0.1;# Parameters of pleasure hedonic adaptation
gammaA=0.05;# Parameters of afflictive affects hedonic adaptation
Map=np.ones(hours)*0;# Approach hedonic motivation
Mev=np.ones(hours)*0;# Avoidance hedonic motivation
############################## Behavioral parameters ###################################
r=np.ones(hours);# probability of success of the approach or the avoidance
#r=np.random.random_sample(hours)# random success
#r[120:167]=0;# two days of failure
q=5; thrP=0.25; thrB=0.5;thrA=0.25/5;# parameters of the motivation laws
gammaPMap=1;gammaPMev=1;#weights of the influence of pleasure et afflictive affects in the motivation
coefexp=0.5;#coefficient that reduce the effect of expectations in case of non contact with stimuli
pcontact=0;#probability of contact

############ Creating one random potential activity per day
activities=np.ones(hours)*0.5*-1;#*(np.exp(-gammaP));#
#activities=np.random.random_sample(hours)*2-1#values of activities for each day between -1 and +1
activitiesdone=np.zeros(hours)
####################### Initialization of the dynamical loop ###########################
criteret=0;compttemps=0;
Ps[0]=0;As[0]=0;#initialization of stimuli pleasure
P[0]=0;A[0]=0;#initialization of the total pleasure
MotivP=1-(P[0]+1)**q/((P[0]+1)**q+(1+thrP)**q);#initilization of the pleasure-based motivation
MotivAA=A[0]**q/(A[0]**q+thrA**q);#initialization of the afflictive affects-based motivation
Map[0]=gammaPMap*MotivP+(1-gammaPMap)*MotivAA;#initialization of the approach motivation
Mev[0]=gammaPMev*MotivP+(1-gammaPMev)*MotivAA;#initialization of the avoidance motivation

########################################################################################
############################ Start of the dynamical loop ###############################
########################################################################################
h=0;d=0;# initialization of time (day and hour)
while compttemps<(hours-1):
  compttemps=compttemps+1;#print(compttemps)
  h=h+1;#incrementation of hour
  ########################### Calculation of the motivation
  Pm=P[compttemps-1]*np.exp(-gammaP);#accumulated pleasure
  Am=A[compttemps-1]*np.exp(-gammaA);#accumulated afflictive affects
  MotivP=1-(Pm+1)**q/((Pm+1)**q+(1+thrP)**q);#calculation of the pleasure-based motivation
  if Pm<-1:#if the accumulated pleasure is lower than 1, motivation is bounded to 1
      MotivP=1;
  MotivAA=Am**q/(Am**q+thrA**q);#calculation of the afflictive-based motivation
  Map[compttemps]=gammaPMap*MotivP+(1-gammaPMap)*MotivAA;#approach motivation
  Mev[compttemps]=gammaPMev*MotivP+(1-gammaPMev)*MotivAA;#avoidance motivation
  if h<17:#individual sleeps 8 hours a day
    ############################# Case of positive stimuli valence #####################
    if activities[compttemps]>0:# we check that stimuli is positive
        expectation=activities[compttemps]+activities[compttemps]*(Map[compttemps]-0.5)*(1-activities[compttemps]);# calculation of the expectation
        if Map[compttemps]>thrB:# case of behavioral activation
            if np.random.random()<r[compttemps-1]:# Success of the approach
                Ps[compttemps]=activities[compttemps]+(activities[compttemps]-expectation);# pleasure caused by the stimuli in the case of approach success
            else:# Failure of the approach
                As[compttemps]=coefexp*expectation;#afflictive affects caused by the failure of the approach
                Ps[compttemps]=-coefexp*expectation;#displeasure cause by the failure of the approach
            activitiesdone[compttemps]=1;
        else: # case of no behavioral activation
            if np.random.random()<pcontact:#contact with the stimuli
                As[compttemps]=0;#no afflictive affects caused by stimuli contact
                Ps[compttemps]=activities[compttemps];#pleasure caused by stimuli contact
            else:#no contact with the stimuli
                As[compttemps]=0;#no afflictive affects (nothing happens)
                Ps[compttemps]=0;#no pleasure (nothing happens)
            activitiesdone[compttemps]=-1;

    ############################# Case of negative stimuli valence #####################
    if activities[compttemps]<0:# we check that stimuli is negative
        expectation=activities[compttemps]+activities[compttemps]*(Mev[compttemps]-0.5)*(1+activities[compttemps]);# calculation of the expectation
        if Mev[compttemps]>thrB:# case of behavioral activation
            if np.random.random()<r[compttemps-1]:# Success of the avoidance
                Ps[compttemps]=-coefexp*expectation;# pleasure caused by the stimuli in the case of avoidance success
            else:# Failure of the avoidance
                As[compttemps]=-activities[compttemps]-(activities[compttemps]-expectation);#afflictive affects caused by the failure of the avoidance
                Ps[compttemps]=activities[compttemps]+(activities[compttemps]-expectation);#displeasure cause by the failure of the approach
        else:
            if np.random.random()<pcontact:#contact with the stimuli
                As[compttemps]=0;#no afflictive affects caused by stimuli contact
                Ps[compttemps]=activities[compttemps];#displeasure caused by stimuli contact
            else:#no contact with the stimuli
                As[compttemps]=0;#no afflictive affects (nothing happens)
                Ps[compttemps]=0;#no pleasure (nothing happens)
  P[compttemps]=Pm+Ps[compttemps];#Calculation of the total pleasure
  A[compttemps]=Am+As[compttemps];#Calculation of the total pleasure
  if h==24:#end of the day
     h=0;#reinitialization of hour
     d=d+1;#day increment

#######################################################################################
###################################### Post-processing ################################
#######################################################################################
# Initialisation of the figure
fig, ax = pl.subplots(2,2)
fig.tight_layout(w_pad=3.0)

# Plot of pleasure and afflictive affects
ax[0,0].plot(P,'b')
ax[0,0].set_xlabel('Time (hour)')
ax[0,0].tick_params(axis='y', colors='blue')
ax[0,0].set_ylabel('Pleasure',color ='blue')
#ax[0,0].set_ylim(-0.3, 0.6)
axt = ax[0,0].twinx()
axt.plot(A,'r')
axt.tick_params(axis='y', colors='red')
axt.set_ylabel('Afflictive affects',color ='red')

# Plot of the avoidance motivation
ax[1,0].plot(Mev,'r')
ax[1,0].set_xlabel('Time (hour)')
ax[1,0].tick_params(axis='y', colors='red')
ax[1,0].set_ylabel('Avoidance motivation',color ='red')
ax[1,0].set_ylim(0, 1)
# Plot of the approach motivation
ax[1,1].plot(Map,'b')
ax[1,1].set_xlabel('Time (hour)')
ax[1,1].tick_params(axis='y', colors='blue')
ax[1,1].set_ylabel('Approach motivation',color ='blue')
ax[1,1].set_ylim(0, 1)

# Plot of happiness
ax[0,1].plot(P-A);
ax[0,1].set_xlabel('Time (hour) ')
ax[0,1].set_ylabel('Happiness')

pl.show()

