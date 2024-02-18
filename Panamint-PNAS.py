# -*- coding: utf-8 -*-
"""
This is for the Rainbow Canyon evolution paper.

@author: GPwgqzhangli
"""

import numpy as np                  
import matplotlib.pyplot as plt
import time
###############################################################################
g=9.81                     ## Acceleration of gravity
Cz=10                      ## Chezy resistance coefficient, Cf=1/Cz**2
Cf=1/Cz**2
nt=1.5                     ## sediment transport formula
fis=1                      ## sediment transport formula
lamda=0.35                 ## alluvial porosity
R=1.65
D=60*(10)**(-3)            ## grain size, [m]
beta=0.015*(10)**(-3)      ## wear coefficient, [m-1] 
up_rate_0=3*(10)**(-3)     ## uplift rate [mm/year→m/yr]
Lmr=1                      ## macro-roughness layer thickness [m]
IF=0.005                   ## flood intermittency
pl=0.05                    ## characteristically low probability on the hypsometric curve
ph=1-pl                    ## characteristically high probability on the hypsometric curve
peta=(1-pl)/(ph-pl)
T_critical=0.0495          ## critical shear number
Bb_value=20                ## channel width at upstream end, m                    
Qw_up=1.45                 ## upstream boundary of water discharge per width, m^2/s  
fb=0.5                     ## side wall feed

'initial value Qw_up, etaa,'
etaa_initial_value=0.5     ## the initial alluvial depth etaa covered the bedrock etab.     
etab_down=0             
L=8000
'''
Slope series and change_value series
'''
Slope_series=np.array([0.0156,0.094])
change_value_series=np.array([7400])
sk_initial=change_value_series[0]      ## the initial knickpoint
slope_qac=0.0156
N_feed=1       
NN=Slope_series.size             
Ss=0.7                                 ## sidewall slope
ETA_0=etaa_initial_value+Slope_series[-1]*(L-change_value_series[-1])+Slope_series[0]*(change_value_series[-1])
H_T_AB=5
###############################################################################
'''
Input Parameters
'''
M=200
dx_bar=1/M
x_all=np.linspace(0,1,2*M+1)
XX=x_all
xx_side=x_all[np.arange(0,2*M+1,2)]      ## To vector
eta_x_initial=x_all[np.arange(1,2*M+1,2)]          
eta_x_initial=np.append(0,eta_x_initial)
eta_x_initial=np.append(eta_x_initial,1) ## To scalar

sk0=(sk_initial-L/2/M)*(2*M)/(2*M-1) 
sk0_initial=sk0
XX_end=L-dx_bar/2*(L-sk0_initial)
L_end_initial=L

dtw=0           ## dt for water part part, 's'
Ntoprint=0      
dta=10**(-4)    ## dt for alluvial part,'yr'
dtb=10**(1)     ## dt for bedrock for water,'yr'
Nprint=200      
NNprint=60     

pprintwant=30 
dwant=6       
Want_to_out=12
show_plot_unit=2   ## for the time scale in plot, 0=second,1=day,2=year
show_unit=2        ## for the time scale in 'cvs'. 0=second,1=day,2=year
##############################################################################
'''
Normal depth H and Fr series
'''
if Slope_series[-1]==0:
    H_normal_series=(Qw_up**2*Cf/np.delete(Slope_series,-1)/g)**(1/3)
else:
    H_normal_series=(Qw_up**2*Cf/Slope_series/g)**(1/3)
    
Fr_normal_series=(Qw_up/H_normal_series)/(g*H_normal_series)**(0.5)
###############################################################################
'''
Qac series at initial state.
'''
if Slope_series[-1]==0:
    T_star_up0_series=((Cz)**2*g)**(-1/3)*(Qw_up)**(2/3)*(np.delete(Slope_series,-1))**(2/3)/R/D
else:
    T_star_up0_series=((Cz)**2*g)**(-1/3)*(Qw_up)**(2/3)*(Slope_series)**(2/3)/R/D
T_star_up00_series=fis*T_star_up0_series-T_critical
T_star_up00_series[np.where(T_star_up00_series<0)[0]]=0          
Qac_up_series=4*D*(R*D*g)**(1/2)*(T_star_up00_series)**nt
'''Qac at upstream end'''
Qac_up=Qac_up_series[0] 
###############################################################################
'''sediment feed''' ## feed sediment at upstream, m2/s              
T_star_slope_qac=((Cz)**2*g)**(-1/3)*(Qw_up)**(2/3)*(slope_qac)**(2/3)/R/D
Qac_slope_qac=4*D*(R*D*g)**(1/2)*(fis*T_star_slope_qac-T_critical)**nt  
Qs_feed=N_feed*Qac_slope_qac
Qs_feed_str=str(round(Qs_feed,6))
###############################################################################
TotalTime=dtb*(365.25*24*60*60)*Nprint*NNprint    ## unit 's'.
everyPlotTime=dtb*(365.25*24*60*60)*Nprint*dwant
everySaveTime=dtb*(365.25*24*60*60)*Nprint*Want_to_out
plot_infor=str(round(TotalTime/(365.25*24*60*60),2))+'yr/'+str(round(everyPlotTime/(365.25*24*60*60),4))+'yr'
save_infor=str(round(TotalTime/(365.25*24*60*60),2))+'yr/'+str(everySaveTime/(365.25*24*60*60))+'yr'
print(plot_infor +'\n'+save_infor)
############x_label#########
dtw_d=str(round(dtw,7))+'s('+str(Ntoprint)+')'
dta_d=str(round(dta,7))+'yrs'
dtb_d=str(round(dtb,7))+'yrs'
labxy_xdxdt='M='+str(M)+'m dtw='+dtw_d+' dta='+dta_d+' dtb='+dtb_d ## for plot of x_label.
#### time part #####
time_s=TotalTime
time_day=time_s/24/60/60
time_year=time_s/365.25/24/60/60
print('time='+str(time_s)+'s='+str(round(time_day,2))+'day='+str(round(time_year,2))+'yr')
time_str='time='+str(round(time_year,2))+'yr'
                           #### time part for 'cvs'#####
time_plot_show=[[]]*3
time_plot_show[0]=str(round(time_s,2))+'s'
time_plot_show[1]=str(round(time_day,2))+'day'
time_plot_show[2]=str(round(time_year,2))+'yrs'
                           #### time part for 'cvs'#####
time_series_total=dtb*(365.25*24*60*60)*Nprint*(NNprint+1)
tima_series_pieces=dtb*(365.25*24*60*60)*Nprint
tima_want_series_pieces=dtb*(365.25*24*60*60)*Nprint*Want_to_out
time_total=np.append(0,np.arange(0,time_series_total,tima_series_pieces))
timewant_s=np.append(0,np.arange(0,time_series_total,tima_want_series_pieces))
timewant_day=timewant_s/24/60/60
timewant_year=timewant_s/365.25/24/60/60
timewant_show=[[]]*3
timewant_show[0]=timewant_s
timewant_show[1]=timewant_day
timewant_show[2]=timewant_year
###############################################################################
HH1=np.zeros([M,M+1])
for i in range(0,M,1):
    HH1[i,i]=0.5
    HH1[i,i+1]=0.5            

HH=np.zeros([M+1,M+2])
for i in range(1,M+1,1):
    HH[i,i]=0.5
    HH[i,i+1]=0.5
HH[0,0]=1

Slope_Eta=np.zeros([M+2,M+2])  
Slope_Eta[0,0]=2
Slope_Eta[0,1]=-2
Slope_Eta[-1,-1]=-2
Slope_Eta[-1,-2]=2
Slope_Eta[1,1]=1
Slope_Eta[1,2]=-1
Slope_Eta[-2,-2]=-1
Slope_Eta[-2,-3]=1
for i in range(2,M,1):
    Slope_Eta[i,i-1]=1/2
    Slope_Eta[i,i+1]=-1/2   
##########################
Bb_initial0=Bb_value*np.ones(len(eta_x_initial))
QWW_up=Qw_up*Bb_value
QFEED_up=Qs_feed*Bb_value

H_2=np.zeros([M+2,NNprint+1])
U_2=np.zeros([M+1,NNprint+1])
Qa_2=np.zeros([M+1,NNprint+1])
etab_2=np.zeros([M+2,NNprint+1])
slope_2=np.zeros([M+2,NNprint+1])

pa_2=np.zeros([M+2,NNprint+1])

Erosion_2=np.zeros([M+2,NNprint+1])
Is_2=np.zeros([M+2,NNprint+1])

sk=np.zeros(NNprint+1)          ## the moving point sk(t)
sk[0]=sk_initial
sk00=np.zeros(NNprint+1)      
L_end00=np.zeros(NNprint+1)   
sk00[0]=sk0_initial
L_end00[0]=L_end_initial
skdot=np.zeros(NNprint+1)       ## the speed of moving point d(sk)/d(t).

Bb_bed_2=Bb_value*np.zeros([M+2,NNprint+1])
eta_top_2=np.zeros([M+2,NNprint+1])
Ss_side_2=Ss*np.ones([M+2,NNprint+1]) 
Bs_side_2=np.zeros([M+2,NNprint+1])

Qw_vary=np.zeros(M+1)
H_vary=np.zeros(M+2)
U_vary=np.zeros(M+1)
Qac_vary=np.zeros(M+1)
Qa_vary=np.zeros(M+1)
etab_vary=np.zeros(M+2)
slope_vary=np.zeros(M+2)
pa_vary=np.zeros(M+2)
E_vary=np.zeros(M+2)
Bb_vary=Bb_value*np.ones(M+2)
eta_top_vary=np.zeros(M+2)
Bs_side_vary=np.zeros(M+2)
Ss_side_vary=np.zeros(M+2)
Is_vary=np.zeros(M+2)

xx0=np.zeros([M+1,NNprint+1])  
xx1=np.zeros([M+2,NNprint+1])  
xx0_vary=np.zeros(M+1)
xx1_vary=np.zeros(M+2)
sk_vary=sk[0]
sk0_vary=sk00[0]
skdot_vary=skdot[0]
sk0_vary=sk00[0]
###############################################################################
'Initial Conditions'
'the initial of bed elevation.'
change_point0_series=np.zeros([NN-1])
change_point0_x_series=np.zeros([NN-1])
for i in range(0,NN-1,1):
    change_point0_series_0=np.where(XX>=change_value_series[i]/(L))
    change_point0_series[i]=int(change_point0_series_0[0][0])
    change_point0_x_series_0=np.where(eta_x_initial>=change_value_series[i]/(L))
    change_point0_x_series[i]=int(change_point0_x_series_0[0][0])

etab_initial=(1-eta_x_initial)*(L-sk0)*Slope_series[1]   
etab_initial[0]=etab_initial[1]
etab_initial[-1]=etab_initial[-2]    
'initial of Qw, H, U, Fr'
Qw_initial=Qw_up*np.ones(M+1)
'the initial of H.'
list_H_series=[[]]*NN
list_H_series[NN-1]=H_normal_series[NN-1]*np.ones(M+2-int(change_point0_x_series[NN-2])) 
for i in range(NN-2,0,-1):
    list_H_series[i]=H_normal_series[i]*np.ones(int(change_point0_x_series[i])-int(change_point0_x_series[i-1])) 
list_H_series[0]=H_normal_series[0]*np.ones(int(change_point0_x_series[0]))
for i in range(0,NN-1,1):
    list_H_series[i+1]=np.append(list_H_series[i],list_H_series[i+1])
H_initial=np.array(list_H_series[i+1])
H_initial[0]=H_initial[1]
H_initial[-1]=H_initial[-2]
'initial of U'
U_initial=Qw_initial/np.dot(HH,H_initial)
'the initial of Qac' 
Qac_initial=Qac_up*np.ones(M+1)
slope_initial=np.dot(Slope_Eta,etab_initial)/dx_bar/(L-sk0)
slope_initial[0]=slope_initial[1]                 
slope_initial[-1]=slope_initial[-2]
'the initial of pa'
pa_initial=Qs_feed/Qac_up_series[1]*np.ones(M+2)
Qa_initial=Qac_initial*pa_initial[1:]
Qa_initial[0]=Qs_feed
E_initial=IF*beta*np.append(Qs_feed,Qa_initial)*(1-pa_initial)
'the initial of cross-section'
Bb_initial=Bb_initial0
Ss_side_initial=Ss*np.ones(M+2) 
xx0_initial=L_end_initial-(1-xx_side)*(L-sk0_initial)
xx1_initial=L_end_initial-(1-eta_x_initial)*(L_end_initial-sk0_initial)
eta_top_initial=-Slope_series[0]*(xx1_initial-xx1_initial[1])+(etab_initial[1]+H_T_AB)
sk_initial_id=np.where(sk_initial==xx1_initial)[0][0]
eta_top_initial[sk_initial_id:]=(etab_initial[sk_initial_id]+H_T_AB)*np.ones(M+2-sk_initial_id)
eta_top_initial[0]=eta_top_initial[1]
eta_top_initial[-1]=eta_top_initial[-2]

Bs_side_initial=(eta_top_initial-etab_initial)/Ss_side_initial
Bs_side_initial[0]=Bs_side_initial[1]
Bs_side_initial[-1]=Bs_side_initial[-2]
Is_initial=np.zeros(M+2)
skdot_initial=0
############################################################################
print('sk at '+str(sk_initial)+'m in ['+str(sk0)+','+str(L)+']' )
############################################################################    
title_p='Rainbow Canyon'+'\n'+time_str
xtitle_p='x (m) \n D='+str(D*10**3)+'mm fb='+str(fb)+ ' $q_{af}$='+Qs_feed_str+'$m^2$/s v='+str(up_rate_0*10**3)+'mm/yr'
#################################################################################     
H_2[:,0]=H_initial
U_2[:,0]=U_initial
Qa_2[:,0]=Qa_initial
etab_2[:,0]=etab_initial
slope_2[:,0]=slope_initial
pa_2[:,0]=pa_initial
Bb_bed_2[:,0]=Bb_initial
Bs_side_2[:,0]=Bs_side_initial
eta_top_2[:,0]=eta_top_initial
Is_2[:,0]=Is_initial

Erosion_2[:,0]=E_initial
sk[0]=sk_initial
sk00[0]=sk0_initial
L_end00[0]=L_end_initial
skdot[0]=skdot_initial
xx1[:,0]=xx1_initial 
xx0[:,0]=xx0_initial

Qw_vary=Qw_initial
H_vary=H_initial
U_vary=U_initial

Qac_vary=Qac_initial
etab_vary=etab_initial
slope_vary=slope_initial
pa_vary=pa_initial
Bb_vary=Bb_initial 
eta_top_vary=eta_top_initial
Ss_side_vary=Ss_side_initial
Bs_side_vary=Bs_side_initial
Is_vary=Is_initial
sk_vary=sk_initial
sk0_vary=sk0_initial
L_end_vary=L_end_initial
skdot_vary=skdot_initial
xx1_vary=xx1_initial
xx0_vary=xx0_initial
Bb_vary_half=np.dot(HH,Bb_vary)

start=time.time()
###############################################################################
for kk in range(0,NNprint,1):
    if kk%pprintwant==0:
         print(kk/pprintwant,'Finish!')
    ###############################################################################
    for k in range(0,Nprint,1):
        L_reach=L_end_vary-sk0_vary
        H_vary=(QWW_up**2*Cf/slope_vary/g/Bb_vary**2)**(1/3)
        H_vary[0]=H_vary[1]
        H_vary[-1]=H_vary[-2]
        
        XX00=sk_vary
        YY00=etab_vary[1]+H_T_AB
        eta_top_vary=YY00-0*(xx1_vary-XX00)
        eta_top_vary[0]=eta_top_vary[1]
        eta_top_vary[-1]=eta_top_vary[-2]
                                
        slope_vary=np.dot(Slope_Eta,etab_vary)/dx_bar/(L_reach)
        slope_vary[0]=slope_vary[1]
        slope_vary[-1]=slope_vary[-2]
        
        Bs_side_vary=(eta_top_vary-etab_vary)/Ss          
        Is_vary=2*fb*Bs_side_vary*E_vary
    ###############################################################################
 ###############################################################################
    H_2[:,kk+1]=H_vary
    Qa_2[:,kk+1]=Qa_vary
    etab_2[:,kk+1]=etab_vary
    slope_2[:,kk+1]=slope_vary
    pa_2[:,kk+1]=pa_vary
    Erosion_2[:,kk+1]=E_vary
    Bs_side_2[:,kk+1]=Bs_side_vary
    eta_top_2[:,kk+1]=eta_top_vary
    Is_2[:,kk+1]=Is_vary
    xx0[:,kk+1]=xx0_vary
    xx1[:,kk+1]=xx1_vary     
    sk[kk+1]=sk_vary
    sk00[kk+1]=sk0_vary
    L_end00[kk+1]=L_end_vary
    skdot[kk+1]=skdot_vary
###############################################################################
elapsed=(time.time()-start)
print('Time Used =', elapsed/60/60, 'h')

WaterSurface=etab_2+H_2

Erosion=np.zeros([M+1,NNprint+1])
Erosion=Erosion_2[0:-1,:]
Erosion[:,0]=np.zeros([M+1])
Erosion_m_yr=Erosion*(365.25*24*60*60)

Is_2_m_yr=Is_2*(365.25*24*60*60)
###############################################################################        
TStop=NNprint+int(dwant)
dwant=int(dwant)
###############################################################################
'Save data Qw,U,H,Qac,Qa,Fr,etaa,etab,Slope,EtaaEtab,WaterSurface,p,pa as .xlsx'
Want_N=int(NNprint/Want_to_out)
list_data=[Qa_2,Erosion_m_yr,H_2,etab_2,slope_2,WaterSurface,pa_2,eta_top_2,Ss_side_2,Is_2_m_yr,Bs_side_2]
list_data00=[xx0,xx1,sk,skdot,sk00,L_end00]
name1=(['qa','Erosion_m_yr'])
name2=(['H','etab','Slope','WaterSurface','p','TopFlat','SideWallSlope','SideFeed_m2_yr','CanyonTopWidth'])
v_N=int(len(name1))
s_N=int(len(name2))  ## v_N+s_N=len(list_data)
list_data1=[[]]*v_N
list_data11=[[]]*v_N
list_data2=[[]]*s_N
list_data21=[[]]*s_N
###############################################################################
'Simple Plots'
pp=1  # pp=0, long title; pp=1, short title. 
xpp=1 # xpp=0 for short x-axis title; pp=1,long x-axis title. 

title_print=[[]]*2
xtitle_print=[[]]*2
title_print[0]='Rainbow Canyon'+'\n'+time_str
title_print[1]=title_p
xtitle_print[0]='reach length x (m)'
xtitle_print[1]=xtitle_p
ytitle_print=(['$q_{a}$ ($m^2$/s)','Erosion (m/yr)',
'H (m)','Bed Elevation $\eta$ (m)',
'Slope S',
'$\eta$+H (m)',
'$p$','Canyon Top Elevation $\eta_{T}$ (m)','$S_s$','Side Wall Feed $I_s$ ($m^2$/yr)','Canyon Top Width $B_{T}$ (m)'])

for i in range(0,v_N,1):
    print(i)
    
    fh=plt.figure(figsize=(18,14))
    ax=plt.gca()
    for j in range(0,TStop,dwant):
        plt.plot(list_data00[0][:,j],list_data[i][:,j],lw=2.5)
    
    colormap=plt.cm.brg
    colors=[colormap(idx) for idx in np.linspace(1,0,len(ax.lines))]
    for idx,line in enumerate(ax.lines):
        line.set_color(colors[idx]) 
        
    plt.plot(list_data00[0][:,0],list_data[i][:,0],'-*g',label='0',lw=3)
    plt.plot(list_data00[0][:,-1],list_data[i][:,-1],'--*k',label=plot_infor,lw=3)
    
    plt.xlabel(xtitle_print[xpp],fontsize=24)
    plt.ylabel(ytitle_print[i],fontsize=24)
    plt.title(title_print[pp],fontsize=24)
    ax.tick_params(labelsize=24)
    plt.legend(loc='upper right',fontsize=24)
    plt.grid()
    plt.savefig(name1[i]+time_plot_show[show_plot_unit]+'.png')
    
for i in range(v_N,v_N+s_N,1):
    print(i)
    
    fh=plt.figure(figsize=(18,14))
    ax=plt.gca()
    for j in range(0,TStop,dwant):
        plt.plot(list_data00[1][:,j],list_data[i][:,j],lw=2.5)
    
    colormap=plt.cm.brg
    colors=[colormap(idx) for idx in np.linspace(1,0,len(ax.lines))]
    for idx,line in enumerate(ax.lines):
        line.set_color(colors[idx]) 
        
    plt.plot(list_data00[1][:,0],list_data[i][:,0],'-*g',label='0',lw=3)
    plt.plot(list_data00[1][:,-1],list_data[i][:,-1],'--*k',label=plot_infor,lw=3)    
    
    plt.xlabel(xtitle_print[xpp],fontsize=24)
    plt.ylabel(ytitle_print[i],fontsize=24)
    plt.title(title_print[pp],fontsize=24)
    ax.tick_params(labelsize=24)
    plt.legend(loc='upper right',fontsize=24)
    plt.grid()
    plt.savefig(name2[i-v_N]+time_plot_show[show_plot_unit]+'.png')

for i in range(v_N+s_N-1,v_N+s_N,1):
    print(i)
    
    fh=plt.figure(figsize=(18,14))
    ax=plt.gca()
    for j in range(0,TStop,dwant):
        plt.plot(list_data00[1][1:-1,j],list_data[i][1:-1,j]/2,lw=2.5)
        plt.plot(list_data00[1][1:-1,j],-list_data[i][1:-1,j]/2,lw=2.5)
    colormap=plt.cm.brg
    colors=[colormap(idx) for idx in np.linspace(1,0,len(ax.lines))]
    for idx,line in enumerate(ax.lines):
        line.set_color(colors[idx]) 
        
    plt.plot(list_data00[1][1:-1,0],list_data[i][1:-1,0]/2,'-*g',label='0',lw=3)
    plt.plot(list_data00[1][1:-1,0],-list_data[i][1:-1,0]/2,'-*g',label='0',lw=3)
    plt.plot(list_data00[1][1:-1,-1],list_data[i][1:-1,-1]/2,'--*k',label=plot_infor,lw=3)    
    plt.plot(list_data00[1][1:-1:,-1],-list_data[i][1:-1,-1]/2,'--*k',label=plot_infor,lw=3)    

    plt.xlabel(xtitle_print[xpp],fontsize=24)
    plt.ylabel(ytitle_print[i],fontsize=24)
    plt.title(title_print[pp],fontsize=24)
    ax.tick_params(labelsize=24)
    plt.legend(loc='upper right',fontsize=24)
    plt.grid()
    plt.savefig(name2[i-v_N]+time_plot_show[show_plot_unit]+'.png')

###############################################################################
for i in range(0,v_N,1):
    print(i)
    list_data11[i]=np.zeros([M+2,Want_N+2])
    list_data11[i][0,:]=timewant_show[show_unit]
    list_data11[i][1:(M+2),1:(Want_N+2)]=list_data[i][:,np.arange(0,(NNprint+1),Want_to_out)]
    np.savetxt(name1[i]+time_plot_show[show_plot_unit]+'.csv',list_data11[i],delimiter=',')
for i in range(0,s_N,1):
    print(i)
    if i==(s_N-1):
        print("Finish!!")    
    list_data21[i]=np.zeros([M+3,Want_N+2])
    list_data21[i][0,:]=timewant_show[show_unit]
    list_data21[i][1:(M+3),1:(Want_N+2)]=list_data[i+v_N][:,np.arange(0,(NNprint+1),Want_to_out)]
    np.savetxt(name2[i]+time_plot_show[show_plot_unit]+'.csv',list_data21[i],delimiter=',')

list_data111=np.zeros([M+2,Want_N+2])
list_data111[0,:]=timewant_show[show_unit]
list_data111[1:(M+2),1:(Want_N+2)]=list_data00[0][:,np.arange(0,(NNprint+1),Want_to_out)]
np.savetxt('x_vector'+time_plot_show[show_plot_unit]+'.csv',list_data111,delimiter=',')

list_data112=np.zeros([M+3,Want_N+2])
list_data112[0,:]=timewant_show[show_unit]
list_data112[1:(M+3),1:(Want_N+2)]=list_data00[1][:,np.arange(0,(NNprint+1),Want_to_out)]
np.savetxt('x_scalar'+time_plot_show[show_plot_unit]+'.csv',list_data112,delimiter=',')

list_data01=np.zeros([2,Want_N+2])
list_data01[0,:]=timewant_show[show_unit]
list_data01[:,0]=np.zeros(2)
list_data01[1:2,1:(Want_N+2)]=list_data00[2][np.arange(0,(NNprint+1),Want_to_out)]
np.savetxt('sk'+time_plot_show[show_plot_unit]+'.csv',list_data01,delimiter=',')

list_data02=np.zeros([2,Want_N+2])
list_data02[0,:]=timewant_show[show_unit]
list_data02[:,0]=np.zeros(2)
list_data02[1:2,1:(Want_N+2)]=list_data00[3][np.arange(0,(NNprint+1),Want_to_out)]
np.savetxt('skSpeed'+time_plot_show[show_plot_unit]+'.csv',list_data02,delimiter=',') 

##############################
''' complex plots'''
fh=plt.figure(figsize=(18,14))
ax=plt.gca()

colormap=plt.cm.brg
colors=[colormap(idx) for idx in np.linspace(1,0,len(np.arange(0,TStop,dwant)))]

colormap=plt.cm.Greys_r
colors1=[colormap(idx) for idx in np.linspace(0.6,0,len(np.arange(0,TStop,dwant)))]

for j in range(0,TStop,dwant):
    i_n=np.where(j==np.arange(0,TStop,dwant))[0][0]
    plt.plot(xx1[1:-1,j],etab_2[1:-1,j],lw=2.5,color=colors[i_n])

XX0=0   
for j in range(0,TStop,dwant):
    i_n=np.where(j==np.arange(0,TStop,dwant))[0][0]
    plt.plot([XX0,xx1[1,j]],[(xx1[1,j]-XX0)*Slope_series[0]+etab_2[1,j],etab_2[1,j]],lw=2.5,color=colors[i_n])
##########    eta_top   #########    
for j in range(0,TStop,dwant):
    i_n=np.where(j==np.arange(0,TStop,dwant))[0][0]
    plt.plot(xx1[1:-1,j],eta_top_2[1:-1,j],':',lw=2.5,color=colors1[i_n])

XX0=0   
for j in range(0,TStop,dwant):
    plt.plot([XX0,xx1[1,j]],[(xx1[1,j]-XX0)*Slope_series[0]+eta_top_2[1,j],eta_top_2[1,j]],':',lw=2.5,color=colors1[i_n])

plt.plot(xx1[1:-1,j],etab_2[1:-1,j],'-*k',label='$\eta_{b}$'+plot_infor,lw=3)

plt.xlabel(xtitle_p,fontsize=24)
plt.ylabel('Elevations',fontsize=24)
plt.title('Rainbow Canyon'+'\n'+time_str,fontsize=24)
ax.tick_params(labelsize=24)
plt.legend(loc='lower left',fontsize=24)
plt.grid()
plt.savefig('Watersurface.png')  

##############################################################################
                               ###'''Finish!'''###
