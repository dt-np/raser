#!/usr/bin/env python3

'''
#  Author: Yuhang Tan <tanyuhang@ihep.ac.cn> 
#  Created [2021-05-18 Tues 14:11] 
#  Based on Raser C++ https://github.com/dt-np/raser
'''


from operator import mod
import matplotlib.pyplot as plt
from array import array
import fenics
import numpy as np
import math
import random
import ROOT
import time
import sys

#defien the detector parameter 
class R2dDetector:
    def __init__(self,l_x,l_y):
        self.l_x = l_x #size
        self.l_y = l_y
    def mesh_step(self,lx_step,ly_step):
        self.n_x = int(self.l_x/lx_step) #mesh step
        self.n_y = int(self.l_y/ly_step)
    def set_para(self,doping,voltage,temperature):
        self.d_neff = doping   #dopingX1e12 cm^-3
        self.v_voltage = voltage   #Voltage
        self.temperature = temperature   #Temperature
        self.t_bin = 50e-12
        self.t_end = 6.0e-9
        self.t_start = 0
        self.n_bin = int((self.t_end-self.t_start)/self.t_bin)
        self.positive_cu = ROOT.TH1F("charge+","Positive Current",self.n_bin,0,self.t_end)
        self.negtive_cu = ROOT.TH1F("charge-","Negative Current",self.n_bin,0,self.t_end)
        self.sum_cu = ROOT.TH1F("charge","Current",self.n_bin,0,self.t_end)

#Calculate the weighting potential and electric field
class Fenics_cal:
    #parameter of SiC
    def __init__(self,my_d, model="PIN"):
        self.l_x = my_d.l_x #size
        self.l_y = my_d.l_y
        self.n_x = my_d.n_x #n step
        self.n_y = my_d.n_y
        self.lx_step = self.l_x/self.n_x
        self.ly_step = self.l_y/self.n_y
        self.model = model
        self.p_electric = []
        self.w_p_electric = []
        #perm_sic = 9.76  #SiCPermittivity
        perm_sic = 11.7   #SiPermittivity
        e0 = 1.60217733e-19
        perm0 = 8.854187817e-12   #F/m
        self.f_value = e0*my_d.d_neff*1e6/perm0/perm_sic
        self.tol = 1e-14
        #fenics space        
        mesh = fenics.RectangleMesh(fenics.Point(0, 0), fenics.Point(my_d.l_x, my_d.l_y), int(my_d.n_x), int(my_d.n_y))
        self.V = fenics.FunctionSpace(mesh, 'P', 1)

    def fenics_p_electric(self,my_d):    #get the electric potential
        # Define boundary condition
        u_D = fenics.Expression('x[1] < tol? det_voltage : 0', degree=2,tol=1E-14,det_voltage=my_d.v_voltage)
        def boundary(x, on_boundary):
            return abs(x[1])<self.tol or abs(x[1]-my_d.l_y)<self.tol
        bc = fenics.DirichletBC(self.V, u_D, boundary)
        # # Define variational problem
        u = fenics.TrialFunction(self.V)
        v = fenics.TestFunction(self.V)

        if self.model == "PIN":
            f = fenics.Constant(self.f_value)
        elif self.model == "SiC_LGAD":
            doping_epr = "(6.1942*(1e14)/sqrt(2*3.1415926)/0.13821*exp(-pow(x[1]-0.67,2))/2/0.13821/0.13821 + (1e13))*((1.60217733e-19)*(1e6)/(8.854187817e-12)/(9.76))*1e-12"
            f = fenics.Expression(doping_epr,degree=2)
        else:
            raise NameError(self.model)
        a = fenics.dot(fenics.grad(u), fenics.grad(v))*fenics.dx
        L = f*v*fenics.dx
        # # Compute solution
        u = fenics.Function(self.V)
        fenics.solve(a == L, u, bc)
        self.electric_value = u.compute_vertex_values()   

    def fenics_p_w_electric(self,my_d):  #get the electric weighting potential
        #####Laplace's equation
        u_w_D = fenics.Expression('x[1] < tol? 0 : 1', degree=2,tol=1E-14)
        def boundary_w(x, on_boundary):
            return abs(x[1])<self.tol or abs(x[1]-my_d.l_y)<self.tol
        bc_w = fenics.DirichletBC(self.V, u_w_D, boundary_w)
        # # Define variational problem
        u_w = fenics.TrialFunction(self.V)
        v_w = fenics.TestFunction(self.V)
        f_w = fenics.Constant(0)
        a_w = fenics.dot(fenics.grad(u_w), fenics.grad(v_w))*fenics.dx
        L_w = f_w*v_w*fenics.dx
        # # Compute solution
        u_w = fenics.Function(self.V)
        fenics.solve(a_w == L_w, u_w, bc_w)
        self.w_electric_value = u_w.compute_vertex_values()

    def change_data_form(self,my_d):
        nx = my_d.n_x+1
        ny = my_d.n_y+1
        self.p_w_electric = [ [] for n in range(nx) ]
        self.p_electric = [ [] for n in range(nx) ]
        self.x_position = [ [] for n in range(nx) ]
        self.y_position = [ [] for n in range(nx) ]
        for j in range(ny):
            for i in range(nx):
                self.x_position[i].append(self.lx_step*(i))
                self.y_position[i].append(self.ly_step*(j))
                if (j==0):
                    self.p_w_electric[i].append(0)
                    self.p_electric[i].append(my_d.v_voltage)
                elif(j==ny-1):
                    self.p_w_electric[i].append(1)
                    self.p_electric[i].append(0)
                else:
                    self.p_w_electric[i].append(self.w_electric_value[i+j*nx])
                    self.p_electric[i].append(self.electric_value[i+j*nx])

    def cal_field(self,my_d):
        nx = my_d.n_x
        ny = my_d.n_y
        self.ex_electric = [ [] for n in range(nx) ]
        self.ey_electric = [ [] for n in range(nx) ]
        self.x_f_position = [ [] for n in range(nx) ]
        self.y_f_position = [ [] for n in range(nx) ]
        for j in range(ny):
            for i in range(nx): 
                self.x_f_position[i].append(0.5*self.lx_step*(2*i+1))
                self.y_f_position[i].append(0.5*self.ly_step*(2*j+1))
                self.ex_electric[i].append((self.p_electric[i+1][j]-self.p_electric[i][j])/self.lx_step)
                self.ey_electric[i].append((self.p_electric[i][j+1]-self.p_electric[i][j])/self.ly_step)

    def cal_point_field(self,px_point,py_point,input_value):
        #Interpolation method 
        rex_value=px_point%self.lx_step
        nx_value=int(px_point/self.lx_step)
        rey_value=py_point%self.ly_step
        ny_value=int(py_point/self.ly_step)
        if(rex_value>self.lx_step/2):
            e_v_x1=rex_value-self.lx_step/2
            nx1_v=nx_value
            nx2_v=nx_value+1
        else:
            e_v_x1=rex_value+self.lx_step/2
            e_v_x2=self.lx_step-e_v_x1
            nx1_v=nx_value-1
            nx2_v=nx_value

        if(rey_value>self.ly_step/2):
            e_v_y1=rey_value-self.ly_step/2
            ny1_v=ny_value
            ny2_v=ny_value+1
        else:
            e_v_y1=rey_value+self.ly_step/2
            e_v_y2=self.ly_step-e_v_y1
            ny1_v=ny_value-1
            ny2_v=ny_value

        if (nx_value<=0):
            r_u=0
            nx1_v=nx2_v
        elif (nx_value>=self.n_x-1):
            r_u=0
            nx2_v=nx1_v
        else:
            r_u=e_v_x1/self.lx_step

        if (ny_value<=0):
            r_t=0
            ny1_v=ny2_v
        elif (ny_value>=self.n_y-1):
            r_t=0
            ny2_v=ny1_v
        else:
            r_t=e_v_y1/self.ly_step

        value_11=input_value[nx1_v][ny1_v]
        value_21=input_value[nx2_v][ny1_v]
        value_12=input_value[nx1_v][ny2_v]
        value_22=input_value[nx2_v][ny2_v]
        out_field=0.0
        out_field=(1-r_u)*(1-r_t)*value_11
        out_field+=r_u*(1-r_t)*value_21
        out_field+=r_u*r_t*value_22
        out_field+=(1-r_u)*r_t*value_12
        return out_field

#####The tracks of the different incident paticle (laser)
class TCTTracks:
    def t_mip(self, i_z, i_r, zlen, rlen, nocarrier):
        self.p_tracks = [ [ [] for n in range(i_z) ] for m in range (i_r)]
        for i in range(i_r):
            p_entry = [0, i]
            p_exit  = [zlen, i]
            p_x=p_entry[0]
            p_y=p_entry[1]
            for j in range(i_z):
                x_div_point = j*(zlen/i_z)+ (zlen/(i_z*2))
                y_div_point = i*(rlen/i_r)+ (rlen/(i_r*2))
                self.p_tracks[i][j].append(x_div_point)
                self.p_tracks[i][j].append(y_div_point)

####The tracks of the different incident paticle (mip)
class Tracks:
    def t_mip(self,p_entry,p_exit,n_div):
        self.p_entry = p_entry
        self.p_exit = p_exit
        self.p_tracks = [ [] for n in range(n_div-1) ]
        self.d_tracks = []
        p_x=p_entry[0]
        p_y=p_entry[1]
        # self.py_entry=[]
        for i in range(n_div-1):
            x_div_point = (p_exit[0]-p_entry[0])/n_div*i+p_entry[0]+(p_exit[0]-p_entry[0])/(2*n_div)
            y_div_point = (p_exit[1]-p_entry[1])/n_div*i+p_entry[1]+(p_exit[1]-p_entry[1])/(2*n_div)
            self.p_tracks[i].append(x_div_point)
            self.p_tracks[i].append(y_div_point)
            self.d_tracks.append(math.sqrt(math.pow(x_div_point-p_x,2)+math.pow(y_div_point-p_y,2)))
            p_x = x_div_point
            p_y = y_div_point
            
# mobility model
def si_mobility(charge,aver_e,my_detector):
    T=my_detector.temperature
    E=aver_e
    Neff=abs(my_detector.d_neff)
    alpha = 0.72*math.pow(T/300.0,0.065)
    if(charge>0):
        ulp = 460.0 * math.pow(T / 300.0, -2.18)
        uminp = 45.0*math.pow(T / 300.0, -0.45)
        Crefp = 2.23e17*math.pow(T / 300.0, 3.2)
        betap = 1.0
        vsatp = 9.05e6 * math.sqrt(math.tanh(312.0/T))
        lfm = uminp + (ulp-uminp)/(1.0 + math.pow(Neff*1e12 / Crefp, alpha))
        hfm = 2*lfm / (1.0+math.pow(1.0 + math.pow(2*lfm * E / vsatp, betap), 1.0 / betap))                       
    else:
        uln = 1430.0 * math.pow(T / 300.0, -2.0)
        uminn = 80.0*math.pow(T / 300.0, -0.45)
        Crefn = 1.12e17*math.pow(T/300.0,3.2)
        betan = 2
        vsatn = 1.45e7 * math.sqrt(math.tanh(155.0/T))
        lfm = uminn + (uln-uminn)/ (1.0 + math.pow(Neff*1e12 / Crefn, alpha))
        hfm = 2*lfm / (1.0+math.pow(1.0 + math.pow(2*lfm * E / vsatn, betan), 1.0/betan))
    return hfm

def sic_mobility(charge,aver_e,my_detector):
    T=my_detector.temperature
    E=aver_e
    Neff=abs(my_detector.d_neff)
    if(charge>0):
        alpha = 0.34
        ulp = 124 * math.pow(T / 300, -2)
        uminp = 15.9
        Crefp = 1.76e19
        betap = 1.213 * math.pow(T / 300.0, 0.17)
        vsatp = 2e7 * math.pow(T / 300.0, 0.52)
        lfm = uminp + ulp/(1.0 + math.pow(Neff*1e12 / Crefp, alpha))
        hfm = lfm / (math.pow(1.0 + math.pow(lfm * E / vsatp, betap), 1.0 / betap))                        
    if(charge<0):
        alpha = 0.61
        ulp = 947 * math.pow(T / 300, -2)
        Crefp = 1.94e19
        betap = 1 * math.pow(T / 300, 0.66)
        vsatp = 2e7 * math.pow(T / 300, 0.87)
        lfm = ulp/ (1 + math.pow(Neff*1e12 / Crefp, alpha))
        hfm = lfm / (math.pow(1.0 + math.pow(lfm * E / vsatp, betap), 1.0/betap))                      
    return hfm 

class Drifts:
    def __init__(self,my_track):
        self.muhh=1650   #mobility related with the magnetic field (now silicon useless)
        self.muhe=310
        self.BB=np.array([0,0])
        self.sstep=1 #drift step
        self.kboltz=8.617385e-5 #eV/K
        self.max_drift_len=1e9 #maximum diftlenght [um]
        self.d_dic_n = {}
        self.d_dic_p = {}
        self.s_time = 0
        for n in range(len(my_track.p_tracks)):
            self.d_dic_n["tk_"+str(n+1)] = [ [] for n in range(4) ]
            self.d_dic_p["tk_"+str(n+1)] = [ [] for n in range(4) ]
    def initial_parameter(self):
        self.end_cond=0
        self.d_time=0
        self.path_len=0
        self.n_step=0
        self.charge=0

    def delta_p(self):
        # magnetic field effect
        if(self.charg)>0:
            FF=self.e_field+self.muhh*np.cross(self.e_field,self.BB)
        else:
            FF=self.e_field-self.muhe*np.cross(self.e_field,self.BB)
        #the delta x with electric field
        if(np.linalg.norm(FF)!=0):
            self.delta_x=-self.sstep*self.charg*FF[0]/np.linalg.norm(FF)
            self.delta_y=-self.sstep*self.charg*FF[1]/np.linalg.norm(FF)
        else:
            self.delta_x=0
            self.delta_y=0

    def drift_v(self,my_d,my_f):
        ex_delta_f = my_f.cal_point_field(self.d_x+self.delta_x,self.d_y+self.delta_y,my_f.ex_electric)
        ey_delta_f = my_f.cal_point_field(self.d_x+self.delta_x,self.d_y+self.delta_y,my_f.ey_electric)
        e_delta_f = np.array([ex_delta_f,ey_delta_f])
        aver_e=(np.linalg.norm(self.e_field)+np.linalg.norm(e_delta_f))/2*1e4
        self.v_drift=sic_mobility(self.charg,aver_e,my_d)*aver_e
        #drift part
        if(self.v_drift==0):
            self.delta_x=0.0
            self.delta_y=0.0
            self.dif_x=0.0
            self.dif_y=0.0
            self.end_cond=9
        else:
            #off when the field gets large enough
            DiffOffField=8  # the silicon value?              
            if(np.linalg.norm(e_delta_f)<DiffOffField):
                self.s_time=self.sstep*1e-4/self.v_drift
                s_sigma=math.sqrt(2*self.kboltz*sic_mobility(self.charg,aver_e,my_d)*my_d.temperature*self.s_time)
                self.dif_x=random.gauss(0,s_sigma)*1e4
                self.dif_y=random.gauss(0,s_sigma)*1e4          
            else:
                self.dif_x=0.0
                self.dif_y=0.0

    def drift_s_step(self,my_d):
        # x axis   
        if((self.d_x+self.delta_x+self.dif_x)>=my_d.l_x): 
            self.d_cx = my_d.l_x
        elif((self.d_x+self.delta_x+self.dif_x)<0):
            self.d_cx = 0
        else:
            self.d_cx = self.d_x+self.delta_x+self.dif_x
        # y axis
        if((self.d_y+self.delta_y+self.dif_y)>=my_d.l_y): 
            self.d_cy = my_d.l_y
        elif((self.d_y+self.delta_y+self.dif_y)<0):
            self.d_cy = 0
        else:
            self.d_cy = self.d_y+self.delta_y+self.dif_y

    def drift_end_condition(self):    
        if(self.wpot>(1-1e-5)):
            self.end_cond=1
        if(self.d_x<=0):
            self.end_cond=2
        if(self.d_y<=0):
            self.end_cond=4
        if(self.path_len>self.max_drift_len):
            self.end_cond=6
        if(self.n_step>10000):
            self.end_cond=7

    def save_inf_track(self,my_d):
        e0 = 1.60217733e-19
        if((self.charge<0 and my_d.v_voltage<0) or (self.charge>0 and my_d.v_voltage>0)): 
            if(self.charg>0):
                self.d_dic_p["tk_"+str(self.n_track)][0].append(self.d_x)
                self.d_dic_p["tk_"+str(self.n_track)][1].append(self.d_y)
                self.d_dic_p["tk_"+str(self.n_track)][2].append(self.charge)
                self.d_dic_p["tk_"+str(self.n_track)][3].append(self.d_time)
            else:
                self.d_dic_n["tk_"+str(self.n_track)][0].append(self.d_x)
                self.d_dic_n["tk_"+str(self.n_track)][1].append(self.d_y)
                self.d_dic_n["tk_"+str(self.n_track)][2].append(self.charge)
                self.d_dic_n["tk_"+str(self.n_track)][3].append(self.d_time)

    def cal_current(self,my_d,my_t):
        my_d.positive_cu.Reset()
        my_d.negtive_cu.Reset()
        my_d.sum_cu.Reset()
        test_p = ROOT.TH1F("test+","test+",my_d.n_bin,0,my_d.t_end)
        test_n = ROOT.TH1F("test-","test-",my_d.n_bin,0,my_d.t_end)
        test_sum = ROOT.TH1F("test sum","test sum",my_d.n_bin,0,my_d.t_end)
        total_pairs=0
        for j in range(len(my_t.p_tracks)):
            for i in range(len(self.d_dic_p["tk_"+str(j+1)][2])):
                test_p.Fill(self.d_dic_p["tk_"+str(j+1)][3][i],self.d_dic_p["tk_"+str(j+1)][2][i])
            test_p = Drifts.get_current_his(self,test_p)           
            for i in range(len(self.d_dic_n["tk_"+str(j+1)][2])):
                test_n.Fill(self.d_dic_n["tk_"+str(j+1)][3][i],self.d_dic_n["tk_"+str(j+1)][2][i])
            test_n = Drifts.get_current_his(self,test_n)
            n_pairs=Drifts.get_energy_loss(self,my_t.d_tracks[j])
            total_pairs+=n_pairs
            test_p.Scale(n_pairs)
            test_n.Scale(n_pairs)            
            my_d.positive_cu.Add(test_p)
            my_d.negtive_cu.Add(test_n)
            test_p.Reset()
            test_n.Reset()
        laudau_t_pairs = Drifts.get_laudau_dis(self,my_t)
        n_scale = laudau_t_pairs/total_pairs
        my_d.positive_cu.Scale(n_scale)
        my_d.negtive_cu.Scale(n_scale)
        my_d.sum_cu.Add(my_d.positive_cu)
        my_d.sum_cu.Add(my_d.negtive_cu)

    def get_current_his(self,histo):
        e0 = 1.60217733e-19
        hist = ROOT.TH1F()
        histo.Copy(hist)
        for i in range(hist.GetNbinsX()):
            histo.SetBinContent(i, hist.GetBinContent(i) \
                /((hist.GetXaxis().GetXmax() - hist.GetXaxis().GetXmin()) \
                / hist.GetNbinsX())*e0)
        return histo

    def get_energy_loss(self,d_drift):
        #Drift distance: d_drift
        loss_energy = 8.4 #silicon carbide /um
        Myf = ROOT.TFile("SiC_5um.root")
        hist = Myf.Get("Edep_device")
        gRandom = ROOT.TRandom3(0)
        ran_energy = hist.GetRandom(gRandom)
        n_pairs = d_drift*ran_energy*1e6/(5*loss_energy)
        return n_pairs

    def get_laudau_dis(self,my_t):
        loss_energy = 8.4
        d_x = math.pow(my_t.p_exit[0]-my_t.p_entry[0],2)
        d_y = math.pow(my_t.p_exit[1]-my_t.p_entry[1],2)
        d_dis = math.sqrt(d_x+d_y)
        d_mpv = 0.04 * math.log(d_dis) + 0.27
        d_FWHM = 0.31 * math.pow(d_dis, -0.17)
        d_da = d_FWHM / 4.
        gRandom = ROOT.TRandom3(0)
        LanConst = gRandom.Landau(d_mpv, d_da)
        if (LanConst > 5.0 * d_mpv):
            LanConst = gRandom.Landau(d_mpv, d_da)
        Lan_pairs = LanConst*1000/loss_energy*d_dis
        return Lan_pairs

    def ionized_drift(self,my_t,my_f,my_d):       
        for i in range(len(my_t.p_tracks)):
            self.n_track=i+1
            for j in range(2):
                if (j==0):
                    self.charg=1 #hole
                if (j==1):
                    self.charg=-1 #electron
                Drifts.initial_parameter(self)
                #generated particles positions
                self.d_x=my_t.p_tracks[i][0]#initial position
                self.d_y=my_t.p_tracks[i][1]
                self.s_time=0.0
                while (self.end_cond==0):
                    if (self.d_y>=(my_d.l_y-1) or self.d_x>=(my_d.l_x-1)):
                        self.end_cond=3  
                    else:                                       
                        # field of the position
                        ex_field = my_f.cal_point_field(self.d_x,self.d_y,my_f.ex_electric)
                        ey_field = my_f.cal_point_field(self.d_x,self.d_y,my_f.ey_electric)
                        self.e_field = np.array([ex_field,ey_field])
                        #delta_poisiton
                        Drifts.delta_p(self)
                        #drift_position
                        Drifts.drift_v(self,my_d,my_f)
                        #drift_next_posiiton
                        Drifts.drift_s_step(self,my_d)
                        #charge_collection
                        delta_Ew=my_f.cal_point_field(self.d_cx,self.d_cy,my_f.p_w_electric)-my_f.cal_point_field(self.d_x,self.d_y,my_f.p_w_electric)
                        self.charge=self.charg*delta_Ew
                        if(self.v_drift!=0):
                            self.d_time=self.d_time+self.s_time
                            self.path_len+=self.sstep
                        self.d_x=self.d_cx
                        self.d_y=self.d_cy
                        self.wpot=my_f.cal_point_field(self.d_x,self.d_y,my_f.p_w_electric)
                        Drifts.save_inf_track(self,my_d)  
                        Drifts.drift_end_condition(self)                                                                         
                    self.n_step+=1
        Drifts.cal_current(self,my_d,my_t)

    def draw_drift_path(self):
        # ROOT.gStyle.SetOptStat(0)
        c1 = ROOT.TCanvas("c1", "canvas1", 200,10,1000, 1000)
        mg = ROOT.TMultiGraph("mg","")
        x_array=array('f')
        y_array=array('f')
        for i in range(len(self.d_dic_p)):
            n=len(self.d_dic_p["tk_"+str(i+1)][0])
            if(n>0):
                x_array.extend(self.d_dic_p["tk_"+str(i+1)][0])
                y_array.extend(self.d_dic_p["tk_"+str(i+1)][1])             
                gr_p = ROOT.TGraph(n,x_array,y_array)
                gr_p.SetMarkerColor(4)
                gr_p.SetLineColor(4)
                gr_p.SetLineStyle(1)
                mg.Add(gr_p)
                del x_array[:]
                del y_array[:]
        for j in range(len(self.d_dic_n)):
            m=len(self.d_dic_n["tk_"+str(j+1)][0])
            if(m>0):
                x_array.extend(self.d_dic_n["tk_"+str(j+1)][0])
                y_array.extend(self.d_dic_n["tk_"+str(j+1)][1])             
                gr_n = ROOT.TGraph(m,x_array,y_array)
                gr_n.SetMarkerColor(2)
                gr_n.SetLineColor(2)
                gr_n.SetLineStyle(1)
                mg.Add(gr_n)
                del x_array[:]
                del y_array[:]
        mg.Draw("APL")
        c1.SaveAs("drift_path.pdf")

#The drift of TCT generated particles
class TCTDrifts:
    def __init__(self,my_track):
        self.muhh=1650   #mobility related with the magnetic field (now silicon useless)
        self.muhe=310
        self.BB=np.array([0,0])
        self.sstep=1 #drift step
        self.kboltz=8.617385e-5 #eV/K
        self.max_drift_len=1e9 #maximum diftlength [um]
        self.d_dic_n = {}
        self.d_dic_p = {}
        self.s_time = 0
        for n in range(len(my_track.p_tracks)*len(my_track.p_tracks[0])):
            self.d_dic_n["tk_"+str(n+1)] = [ [] for n in range(4) ]
            self.d_dic_p["tk_"+str(n+1)] = [ [] for n in range(4) ]

    def initial_parameter(self):
        self.end_cond=0
        self.d_time=0
        self.path_len=0
        self.n_step=0
        self.charge=0

    def delta_p(self):
        # magnetic field effect
        if(self.charg)>0:
            FF=self.e_field+self.muhh*np.cross(self.e_field,self.BB)
        else:
            FF=self.e_field-self.muhe*np.cross(self.e_field,self.BB)
        #the delta x with electric field
        if(np.linalg.norm(FF)!=0):
            self.delta_x=-self.sstep*self.charg*FF[0]/np.linalg.norm(FF)
            self.delta_y=-self.sstep*self.charg*FF[1]/np.linalg.norm(FF)                                                                                                                       
        else:
            self.delta_x=0
            self.delta_y=0

    def drift_v(self,my_d,my_f):
        ex_delta_f = my_f.cal_point_field(self.d_x+self.delta_x, self.d_y+self.delta_y, my_f.ex_electric)
        ey_delta_f = my_f.cal_point_field(self.d_x+self.delta_x, self.d_y+self.delta_y, my_f.ey_electric)
        e_delta_f = np.array([ex_delta_f,ey_delta_f])
        aver_e=(np.linalg.norm(self.e_field)+np.linalg.norm(e_delta_f))/2*1e4
        self.v_drift=si_mobility(self.charg,aver_e,my_d)*aver_e
        #drift part
        if(self.v_drift==0):
            self.delta_x=0.0
            self.delta_y=0.0
            self.dif_x=0.0
            self.dif_y=0.0
            self.end_cond=9
        else:
            #off when the field gets large enough
            DiffOffField=8   # the silicon value ???
            if(np.linalg.norm(e_delta_f)<DiffOffField):
                self.s_time=self.sstep*1e-4/self.v_drift
                s_sigma=math.sqrt(2 *self.kboltz *si_mobility(self.charg, aver_e, my_d) *my_d.temperature *self.s_time)
                self.dif_x=random.gauss(0,s_sigma)*1e4
                self.dif_y=random.gauss(0,s_sigma)*1e4
            else:
                self.dif_x=0.0
                self.dif_y=0.0

    def drift_s_step(self,my_d):
        # x axis   
        if((self.d_x+self.delta_x+self.dif_x)>=my_d.l_x):
            self.d_cx = my_d.l_x
        elif((self.d_x+self.delta_x+self.dif_x)<0):
            self.d_cx = 0
        else:
            self.d_cx = self.d_x+self.delta_x+self.dif_x
        # y axis
        if((self.d_y+self.delta_y+self.dif_y)>=my_d.l_y): 
            self.d_cy = my_d.l_y
        elif((self.d_y+self.delta_y+self.dif_y)<0):
            self.d_cy = 0
        else:
            self.d_cy = self.d_y+self.delta_y+self.dif_y

    def drift_end_condition(self):    
        if(self.wpot>(1-1e-5)):
            self.end_cond=1
        if(self.d_x<=0):
            self.end_cond=2
        if(self.d_y<=0):
            self.end_cond=4
        if(self.path_len>self.max_drift_len):
            self.end_cond=6
        if(self.n_step>10000):
            self.end_cond=7

    def save_inf_track(self,my_d):
        e0 = 1.60217733e-19
        if((self.charge<0 and my_d.v_voltage<0) or (self.charge>0 and my_d.v_voltage>0)):
            if(self.charg>0):
                self.d_dic_p["tk_"+str(self.n_track)][0].append(self.d_x)
                self.d_dic_p["tk_"+str(self.n_track)][1].append(self.d_y)
                self.d_dic_p["tk_"+str(self.n_track)][2].append(self.charge)
                self.d_dic_p["tk_"+str(self.n_track)][3].append(self.d_time)
            else:
                self.d_dic_n["tk_"+str(self.n_track)][0].append(self.d_x)
                self.d_dic_n["tk_"+str(self.n_track)][1].append(self.d_y)
                self.d_dic_n["tk_"+str(self.n_track)][2].append(self.charge)
                self.d_dic_n["tk_"+str(self.n_track)][3].append(self.d_time)

    def cal_current(self,my_d,my_t,nocarrier):
        my_d.positive_cu.Reset()
        my_d.negtive_cu.Reset()
        my_d.sum_cu.Reset()
        test_p = ROOT.TH1F("test+","test+",my_d.n_bin,0,my_d.t_end)
        test_n = ROOT.TH1F("test-","test-",my_d.n_bin,0,my_d.t_end)
        test_sum = ROOT.TH1F("test sum","test sum",my_d.n_bin,0,my_d.t_end)
        for i_r in range(len(my_t.p_tracks)):
            for i_z in range(len(my_t.p_tracks[i_r])):
                for i in range(len(self.d_dic_p["tk_"+str(i_r*(len(my_t.p_tracks[0]))+i_z+1)][2])):
                    test_p.Fill(self.d_dic_p["tk_"+str(i_r*(len(my_t.p_tracks[0]))+i_z+1)][3][i],self.d_dic_p["tk_"+str(i_r*(len(my_t.p_tracks[0]))+i_z+1)][2][i])
                test_p = TCTDrifts.get_current_his(self,test_p)
                test_p.Scale(nocarrier[i_z, i_r])
                my_d.positive_cu.Add(test_p)
                for i in range(len(self.d_dic_n["tk_"+str(i_r*(len(my_t.p_tracks[0]))+i_z+1)][2])):
                    test_n.Fill(self.d_dic_n["tk_"+str(i_r*(len(my_t.p_tracks[0]))+i_z+1)][3][i],self.d_dic_n["tk_"+str(i_r*(len(my_t.p_tracks[0]))+i_z+1)][2][i])
                test_n = TCTDrifts.get_current_his(self,test_n)
                test_n.Scale(nocarrier[i_z, i_r])
                my_d.negtive_cu.Add(test_n)
                test_p.Reset()
                test_n.Reset()
        my_d.sum_cu.Add(my_d.positive_cu)
        my_d.sum_cu.Add(my_d.negtive_cu)
    def get_current_his(self,histo):
        e0 = 1.60217733e-19
        hist = ROOT.TH1F()
        histo.Copy(hist)
        for i in range(hist.GetNbinsX()):
            histo.SetBinContent(i, hist.GetBinContent(i) \
                /((hist.GetXaxis().GetXmax() - hist.GetXaxis().GetXmin()) \
                / hist.GetNbinsX())*e0)
        return histo

    def tct_ionized_drift(self,my_t,my_f,my_d,nocarrier):  
        self.n_track=0
        for n in range(len(my_t.p_tracks)):
            for i in range(len(my_t.p_tracks[n])):
                self.n_track = self.n_track + 1
                for j in range(2):
                    if (j==0):
                        self.charg=1 #hole
                    if (j==1):
                        self.charg=-1 #electron
                    TCTDrifts.initial_parameter(self)
                    #generated particles positions
                    self.d_x=my_t.p_tracks[n][i][0]#initial position
                    self.d_y=my_t.p_tracks[n][i][1]
                    self.s_time=0.0
                    while (self.end_cond==0):
                        if (self.d_y>(my_d.l_y-0.5) or self.d_x>(my_d.l_x-5)):
                            self.end_cond=3  
                        else:                            
                            # field of the position
                            ex_field = my_f.cal_point_field(self.d_x,self.d_y,my_f.ex_electric)
                            ey_field = my_f.cal_point_field(self.d_x,self.d_y,my_f.ey_electric)
                            self.e_field = np.array([ex_field,ey_field])
                            #delta_poisiton
                            TCTDrifts.delta_p(self)
                            #drift_position
                            TCTDrifts.drift_v(self,my_d,my_f)
                            #drift_next_posiiton
                            TCTDrifts.drift_s_step(self,my_d)
                            #charge_collection
                            delta_Ew=my_f.cal_point_field(self.d_cx,self.d_cy,my_f.p_w_electric)-my_f.cal_point_field(self.d_x,self.d_y,my_f.p_w_electric)
                            self.charge=self.charg*delta_Ew
                            if(self.v_drift!=0):
                                self.d_time=self.d_time+self.s_time
                                self.path_len+=self.sstep
                            self.d_x=self.d_cx
                            self.d_y=self.d_cy
                            self.wpot=my_f.cal_point_field(self.d_x,self.d_y,my_f.p_w_electric)
                            TCTDrifts.save_inf_track(self,my_d)  
                            TCTDrifts.drift_end_condition(self)                                                                         
                        self.n_step+=1
        TCTDrifts.cal_current(self,my_d,my_t,nocarrier)

    def draw_drift_path(self):
        c1 = ROOT.TCanvas("c1", "canvas1", 200, 10, 1000, 1000)
        mg = ROOT.TMultiGraph("mg","")
        x_array=array('f')
        y_array=array('f')
        for i in range(len(self.d_dic_p)):
            n=len(self.d_dic_p["tk_"+str(i+1)][0])
            if(n>0):
                x_array.extend(self.d_dic_p["tk_"+str(i+1)][0])
                y_array.extend(self.d_dic_p["tk_"+str(i+1)][1])    
                gr_p = ROOT.TGraph(n,x_array,y_array)
                gr_p.SetMarkerColor(2)
                gr_p.SetLineColor(2)
                gr_p.SetLineStyle(1)
                mg.Add(gr_p)
                del x_array[:]
                del y_array[:]
        for j in range(len(self.d_dic_n)):
            m=len(self.d_dic_n["tk_"+str(j+1)][0])
            if(m>0):
                x_array.extend(self.d_dic_n["tk_"+str(j+1)][0])
                y_array.extend(self.d_dic_n["tk_"+str(j+1)][1])             
                gr_n = ROOT.TGraph(m,x_array,y_array)
                gr_n.SetMarkerColor(4)
                gr_n.SetLineColor(4)
                gr_n.SetLineStyle(1)
                mg.Add(gr_n)
                del x_array[:]
                del y_array[:]
        mg.Draw("APL")
        c1.SaveAs("drift_path.pdf")

class Amplifier():
    def CSA_amp(self,my_d,t_rise,t_fall,trans_imp):
        hist = ROOT.TH1F()
        my_d.sum_cu.Copy(hist)
        max_num = hist.GetNbinsX()
        preamp_Q = [0.0]*max_num
        itot = [0.0]*max_num
        shaper_out_Q = [0.0]*max_num
        shaper_out_V = [0.0]*max_num
        qtot = 0.0
        sh_max = 0.0

        tau_rise = t_rise / 2.2*1e-9
        tau_fall = t_fall / 2.2*1e-9   
        if (tau_rise == tau_fall):
            tau_rise *= 0.9
        time_unit = (hist.GetXaxis().GetXmax()- hist.GetXaxis().GetXmin()) / hist.GetNbinsX()
        for i in range(max_num-1):
            if (i>0):
                preamp_Q[i] = 0.0
                itot[i] = hist.GetBinContent(i)
                preamp_Q[i] = itot[i] * time_unit
                qtot += itot[i] *time_unit
            elif (i==0):
                preamp_Q[i]==0.0
        for i in range(max_num-1):
            dif_shaper_Q = preamp_Q[i]
            if(dif_shaper_Q != 0):
                for j in range(max_num-i):
                    shaper_out_Q[i+j] += tau_fall/(tau_fall+tau_rise)*dif_shaper_Q \
                        *(math.exp(-j*time_unit/tau_fall)-math.exp(-j*time_unit/tau_rise))
            if (abs(shaper_out_Q[i]) > abs(sh_max)):
                sh_max = shaper_out_Q[i]
        for i in range(max_num):
            if sh_max !=0:
                shaper_out_V[i] = shaper_out_Q[i] * trans_imp * 1e15 *qtot / sh_max
                hist.SetBinContent(i,shaper_out_V[i])
            else:
                continue
        return qtot,hist

def draw_plot(my_detector,ele_current,qtot,my_drift):
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas("c", "canvas", 200,10,1000, 1000)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.12)
    # c.SetTopMargin(0.12)
    c.SetBottomMargin(0.14)
    my_detector.sum_cu.GetXaxis().SetTitleOffset(1.2)
    my_detector.sum_cu.GetXaxis().SetTitleSize(0.05)
    my_detector.sum_cu.GetXaxis().SetLabelSize(0.04)
    my_detector.sum_cu.GetXaxis().SetNdivisions(510)
    my_detector.sum_cu.GetYaxis().SetTitleOffset(1.1)
    my_detector.sum_cu.GetYaxis().SetTitleSize(0.05)
    my_detector.sum_cu.GetYaxis().SetLabelSize(0.04)
    my_detector.sum_cu.GetYaxis().SetNdivisions(505)
    my_detector.sum_cu.GetXaxis().CenterTitle()
    my_detector.sum_cu.GetYaxis().CenterTitle()
    my_detector.sum_cu.GetXaxis().SetTitle("Time [s]")
    my_detector.sum_cu.GetYaxis().SetTitle("Current [A]")
    my_detector.sum_cu.Draw("HIST")
    my_detector.positive_cu.Draw("SAME HIST")
    my_detector.negtive_cu.Draw("SAME HIST")
    c.Update()
    if ele_current.GetMinimum() <0:
        rightmax = 1.1*ele_current.GetMinimum()
    else:
        rightmax = 1.1*ele_current.GetMaximum()
    if rightmax == 0:
        n_scale = 0
    elif ele_current.GetMinimum() <0:
        n_scale = ROOT.gPad.GetUymin() /rightmax
    else:
        n_scale = ROOT.gPad.GetUymax() /rightmax
    ele_current.Scale(n_scale)
    ele_current.Draw("SAME HIST")
    my_detector.sum_cu.SetLineColor(3)
    my_detector.positive_cu.SetLineColor(2)
    my_detector.negtive_cu.SetLineColor(4)
    ele_current.SetLineColor(6)
    my_detector.sum_cu.SetLineWidth(2)
    my_detector.positive_cu.SetLineWidth(2)
    my_detector.negtive_cu.SetLineWidth(2)
    ele_current.SetLineWidth(2)    
    axis = ROOT.TGaxis(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(), rightmax, 0, 510, "+L")
    axis.SetLineColor(6)
    axis.SetTextColor(6)
    axis.SetTextSize(0.02)
    axis.SetLabelColor(6)
    axis.SetLabelSize(0.02)
    axis.SetTitle("Ampl [mV]")
    axis.CenterTitle()
    axis.Draw("same")

    legend = ROOT.TLegend(0.5, 0.3, 0.8, 0.6)
    legend.AddEntry(my_detector.negtive_cu, "electron", "l")
    legend.AddEntry(my_detector.positive_cu, "hole", "l")
    legend.AddEntry(my_detector.sum_cu, "e+h", "l")
    legend.AddEntry(ele_current, "electronics", "l")
    legend.SetBorderSize(0)
    legend.SetTextFont(43)
    legend.SetTextSize(40)
    legend.Draw("same")

    c.Update()
    c.SaveAs("basic_infor.pdf")
    
    charge_t=my_detector.sum_cu.Integral() \
        * ((my_detector.sum_cu.GetXaxis().GetXmax() \
        - my_detector.sum_cu.GetXaxis().GetXmin()) \
        / my_detector.sum_cu.GetNbinsX()) * 1e15
    #print(charge_t)
    #print(qtot*1e15)
    my_drift.draw_drift_path()

def TCTdraw_plot(my_detector,ele_current,qtot,my_drift,my_track,nocarrier):
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas("c", "canvas", 200,10,1000, 1000)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.12)
    c.SetTopMargin(0.12)
    c.SetBottomMargin(0.14)
    my_detector.sum_cu.GetXaxis().SetTitleOffset(1.2)
    my_detector.sum_cu.GetXaxis().SetTitleSize(0.05)
    my_detector.sum_cu.GetXaxis().SetLabelSize(0.04)
    my_detector.sum_cu.GetXaxis().SetNdivisions(510)
    my_detector.sum_cu.GetYaxis().SetTitleOffset(1.1)
    my_detector.sum_cu.GetYaxis().SetTitleSize(0.05)
    my_detector.sum_cu.GetYaxis().SetLabelSize(0.04)
    my_detector.sum_cu.GetYaxis().SetNdivisions(505)
    my_detector.sum_cu.GetXaxis().CenterTitle()
    my_detector.sum_cu.GetYaxis().CenterTitle()
    my_detector.sum_cu.GetXaxis().SetTitle("Time [s]")
    my_detector.sum_cu.GetYaxis().SetTitle("Current [A]")
    my_detector.sum_cu.Draw("HIST")
    my_detector.positive_cu.Draw("SAME HIST")
    my_detector.negtive_cu.Draw("SAME HIST")
    c.Update()
    if ele_current.GetMinimum() <0:
        rightmax = 1.1*ele_current.GetMinimum()
    else:
        rightmax = 1.1*ele_current.GetMaximum()
    if rightmax == 0:
        n_scale = 0
    elif ele_current.GetMinimum() <0:
        n_scale = ROOT.gPad.GetUymin() /rightmax
    else:
        n_scale = ROOT.gPad.GetUymax() /rightmax
    ele_current.Scale(n_scale)
    ele_current.Draw("SAME HIST")
    my_detector.sum_cu.SetLineColor(3)
    my_detector.positive_cu.SetLineColor(2)
    my_detector.negtive_cu.SetLineColor(4)
    ele_current.SetLineColor(6)
    my_detector.sum_cu.SetLineWidth(2)
    my_detector.positive_cu.SetLineWidth(2)
    my_detector.negtive_cu.SetLineWidth(2)
    ele_current.SetLineWidth(2)    
    axis = ROOT.TGaxis(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(), rightmax, 0, 510, "+L")
    axis.SetLineColor(6)
    axis.SetTextColor(6)
    axis.SetTextSize(0.02)
    axis.SetLabelColor(6)
    axis.SetLabelSize(0.02)
    axis.SetTitle("Ampl [mV]")
    axis.CenterTitle()
    axis.Draw("same")

    legend = ROOT.TLegend(0.5, 0.3, 0.8, 0.6)
    legend.AddEntry(my_detector.negtive_cu, "electron", "l")
    legend.AddEntry(my_detector.positive_cu, "hole", "l")
    legend.AddEntry(my_detector.sum_cu, "e+h", "l")
    legend.AddEntry(ele_current, "electronics", "l")
    legend.SetBorderSize(0)
    legend.SetTextFont(43)
    legend.SetTextSize(40)
    legend.Draw("same")

    c.Update()
    c.SaveAs("basic_infor.pdf")
    
    charge_t=my_detector.sum_cu.Integral() \
        * ((my_detector.sum_cu.GetXaxis().GetXmax() \
        - my_detector.sum_cu.GetXaxis().GetXmin()) \
        / my_detector.sum_cu.GetNbinsX()) * 1e15
    print(charge_t)
    print(qtot*1e15)
    my_drift.draw_drift_path()

    elecu = [ []for n in range(ele_current.GetNbinsX())]
    for i in range(ele_current.GetNbinsX()):
        elecu[i].append(ele_current.GetBinContent(i))
    mincu = min(elecu)
    return(mincu)

def draw_nocarrier(my_t, my_d, nocarrier):
    c1 = ROOT.TCanvas("c1","canvas2",200,10,1000,1000)
    h = ROOT.TH2D("h","pairs of carrier generation", len(my_t.p_tracks[0]), 0, my_d.l_x, len(my_t.p_tracks), 0, my_d.l_y)
    for i_r in range(len(my_t.p_tracks)):
        for i_z in range(len(my_t.p_tracks[i_r])):
            n = nocarrier[i_z, i_r]
            h.Fill(i_z*(my_d.l_x/len(my_t.p_tracks[0])), i_r*(my_d.l_y/len(my_t.p_tracks)), n)
    h.Draw()
    c1.SaveAs("nocarrier.pdf")

class SPAGeneration():
    def __init__(self):
        ### for Si 1064 IR
        self.tau = 350e-12   #ps
        self.alfa = 987      #SI
        self.power = 1e-11   #J/s
        self.wavelength = 1.064e-6  #m
        self.widthBeamWaist = 5e-6  #m
        self.refractionIndex = 3.51

    def getWidthSquared(self,z):
        return (self.widthBeamWaist**2)*(1+((self.wavelength*z)/(np.pi* (self.widthBeamWaist**2)*self.refractionIndex))**2)
    def getWidth(self,z):
        return np.sqrt(self.getWidthSquared(z))
    def getIntensity(self,r,z,z_o = 0 ):
        widthSquared= self.getWidthSquared(z-z_o)
        intensity = ((2*self.power)/(np.pi*widthSquared))*np.exp((-2*(r**2)/(widthSquared)))*np.exp(-self.alfa*z)
        return intensity
    def getCarrierDensity(self,r,z,z_o = 0 ):
        I = self.getIntensity(r,z,z_o)
        return (self.alfa*I)/(3.6*1.60217657e-19)

class ProjGrid():
    def __init__(self,r_min,r_max,r_nBins,z_max,z_min,z_nBins,z_o):
        self.r_min = r_min
        self.r_max = r_max
        self.r_nBins = r_nBins
        self.z_min = z_min
        self.z_max = z_max
        self.z_nBins = z_nBins
        self.z_o = z_o
        self.rArray = np.linspace(r_min*1e-6,r_max*1e-6,r_nBins)
        self.zArray = np.linspace(z_min*1e-6,z_max*1e-6,z_nBins)
        self.r_step = abs(r_max*1e-6-r_min*1e-6)/r_nBins
        self.z_step = abs(z_max*1e-6-z_min*1e-6)/z_nBins
        self.xArray = self.rArray
        self.yArray = self.rArray
        self.x_step = self.r_step
        self.y_step = self.r_step

class Matplt:
    def plot_basic_info(self,my_f,my_drift):
        plt.figure(figsize=(20,20))

        plt.subplot(2,2,1)
        plt.title('Electric field')
        plt.xlabel('depth [um]')
        plt.ylabel('Electric field [V/um]')
        plt.plot(my_f.y_f_position[26],my_f.ey_electric[26])

        plt.subplot(2,2,2)
        plt.title('Electric field')
        plt.xlabel('X [um]')
        plt.ylabel('Electric field [V/um]')
        plt.plot(my_f.x_f_position[1],my_f.ex_electric[1])

        plt.subplot(2,2,3)
        plt.title('weighting potential')
        plt.xlabel('depth [um]')
        plt.ylabel('Electric potential [V]')
        plt.plot(my_f.y_position[0], my_f.p_w_electric[0])

        plt.subplot(2,2,4)
        plt.title('potential')
        plt.xlabel('depth [um]')
        plt.ylabel('Electric potential [V]')
        plt.plot(my_f.y_position[0], my_f.p_electric[0])
        plt.savefig("test_electric.pdf")
        # plt.show()

class NNocarrier():
    def __init__(self, r, rlen, zlen, rbin, zbin):
        self.z_o = 0   # z focus position in um
        self.rlen = rlen
        self.zlen = zlen
        self.r_min = -r # um
        self.r_max = -r+ rlen # um
        self.r_nBins = rbin
        self.z_min = 0. # um
        self.z_max = zlen # um
        self.z_nBins = zbin

    def Nocarrier(self):
        my_pro = ProjGrid(self.r_min, self.r_max, self.r_nBins, self.z_max, self.z_min, self.z_nBins, self.z_o)
        rGrid, zGrid = np.meshgrid(my_pro.rArray,my_pro.zArray)
        carriergeneration = SPAGeneration()
        CGrid = carriergeneration.getCarrierDensity(rGrid, zGrid, self.z_o*1e-6)
        xGrid = rGrid.copy()
        projGrid = CGrid.copy()

        for i_z in range(self.z_nBins):
            for i_r in range(self.r_nBins):
                x_value = xGrid[i_z, i_r]
                r_valueArray = np.sqrt(x_value*x_value+ my_pro.yArray*my_pro.yArray)
                # Get z value
                z_value = zGrid[i_z, i_r]
                z_valueArray = np.ones_like(r_valueArray)*z_value
                carr_den_projArray = carriergeneration.getCarrierDensity(r_valueArray, z_valueArray, self.z_o*1e-6)
                # Project sum and take into account integral step
                projGrid[i_z, i_r] = (carr_den_projArray.sum()*my_pro.y_step)/2
        projGrid = projGrid*my_pro.x_step*my_pro.z_step
        return(projGrid)


### get the 2D TCT simulation basics information
def twoD_TCT(model="PIN"):
    ### define the structure of the detector
    rlen = 50
    zlen = 1300
    rbin = 50
    zbin = 130
    my_detector = R2dDetector(zlen, rlen)
    my_detector.mesh_step(1,1)
    # my_detector.mesh_step(10,1)
    my_detector.set_para(doping=1.79, voltage=-200, temperature=300)      #doping:e12
    ### get the electric field and weighting potential
    if model == "PIN":
        my_field = Fenics_cal(my_detector)
    
    elif model == "SiC_LGAD":
        my_field = Fenics_cal(my_detector,model="SiC_LGAD")

    else:
        raise NameError(model)

    my_field.fenics_p_electric(my_detector) 
    my_field.fenics_p_w_electric(my_detector)
    my_field.change_data_form(my_detector)
    my_field.cal_field(my_detector)

    num = 70
    cvmin = ROOT.TCanvas("cvmin", "vmin", 100, 100, 800, 600)
    mg = ROOT.TMultiGraph("mg","")
    x_array=array('f')
    y_array=array('f')
    x = [[]for n in range(num+1)]
    y = [ ]
    for i in range(num+1):
        r = -10+ i*(70/num)
        carriergen = NNocarrier(r, rlen, zlen, rbin, zbin)    #Input the position of the laser
        nocarrier = carriergen.Nocarrier()
        ### define the tracks and type of incident particles
        my_track = TCTTracks()
        my_track.t_mip(zbin, rbin, zlen, rlen, nocarrier)
        ### drift of ionized particles
        my_drift = TCTDrifts(my_track)
        my_drift.tct_ionized_drift(my_track,my_field,my_detector,nocarrier)
        ### after the electronics
        my_electronics = Amplifier()
        qtot,ele_current=my_electronics.CSA_amp(my_detector,t_rise=0.4,t_fall=0.2,trans_imp=10)
        ### matlab plot and show
        my_plot = Matplt()
        my_plot.plot_basic_info(my_field,my_drift)
        ### root plot
        draw_nocarrier(my_track, my_detector, nocarrier)
        cumin = TCTdraw_plot(my_detector,ele_current,qtot,my_drift,my_track,nocarrier)
        vmin = [i * 50 for i in cumin]
        y.append(vmin)
        x[i].append(r)
        x_array.extend(x[i])
        y_array.extend(vmin)
        n = 1
        gr = ROOT.TGraph(n,x_array,y_array)
        gr.SetMarkerColor(47)
        gr.SetMarkerStyle(42)
        gr.SetMarkerSize(3)
        mg.Add(gr)
        del x_array[:]
        del y_array[:]
    mg.Draw("APL")
    cvmin.SaveAs("vmin200V.pdf")

    outputfile = ROOT.TFile("dataout200V.root","RECREATE")
    tree = ROOT.TTree("tree","tree")
    v = array('f', [0.])
    z = array('f', [0.])
    tree.Branch("v", v, "v/F")
    tree.Branch("z", z, "z/F")
    for i in range(num+1):
        v[0] = y[i][0]
        z[0] = x[i][0]
        tree.Fill()
    tree.Write()
    outputfile.Close()

    # #Input the position of the laser
    # # define the tracks and type of incident particles
    # carriergen = NNocarrier(25, rlen, zlen, rbin, zbin)    #Input the position of the laser
    # nocarrier = carriergen.Nocarrier() 
    # my_track = TCTTracks()
    # my_track.t_mip(zbin, rbin, zlen, rlen, nocarrier)
    # ## drift of ionized particles
    # my_drift = TCTDrifts(my_track)
    # my_drift.tct_ionized_drift(my_track,my_field,my_detector,nocarrier)
    # ## after the electronics
    # my_electronics = Amplifier()
    # qtot,ele_current=my_electronics.CSA_amp(my_detector,t_rise=0.4,t_fall=0.2,trans_imp=10)
    # ## matlab plot and show
    # my_plot = Matplt()
    # my_plot.plot_basic_info(my_field,my_drift)
    # ## root plot
    # draw_nocarrier(my_track, my_detector, nocarrier)
    # cumin = TCTdraw_plot(my_detector,ele_current,qtot,my_drift,my_track,nocarrier) 
    # print(cumin)

### get the 2D simulation basics information
def twoD_time(model="PIN"):
    ### define the structure of the detector
    my_detector = R2dDetector(100,50)
    my_detector.mesh_step(1,1)
    my_detector.set_para(doping=-10,voltage=-200,temperature=300)
    ### get the electric field and weighting potential
    if model == "PIN":
        my_field = Fenics_cal(my_detector)

    elif model == "SiC_LGAD":
        my_field = Fenics_cal(my_detector,model="SiC_LGAD")

    else:
        raise NameError(model)
    my_field.fenics_p_electric(my_detector) 
    my_field.fenics_p_w_electric(my_detector)
    my_field.change_data_form(my_detector)
    my_field.cal_field(my_detector)  

    ### define the tracks and type of incident particles
    my_track = Tracks()
    my_track.t_mip([50,0],[50,50],100)
    ### drift of ionized particles
    my_drift = Drifts(my_track)
    my_drift.ionized_drift(my_track,my_field,my_detector)
    ### after the electronics
    my_electronics = Amplifier()
    qtot,ele_current=my_electronics.CSA_amp(my_detector,t_rise=0.4,t_fall=0.2,trans_imp=10)
    ### matlab plot and show
    my_plot = Matplt()
    my_plot.plot_basic_info(my_field,my_drift)
    ### root plot
    draw_plot(my_detector,ele_current,qtot,my_drift)

### get the 2D time resolution
def twoD_time_scan(output,number):
    ### define the structure of the detector
    my_detector = R2dDetector(100,100)
    my_detector.mesh_step(25,25)
    my_detector.set_para(doping=-10,voltage=-500,temperature=300)
    ### get the electric field and weighting potential
    my_field = Fenics_cal(my_detector)
    my_field.fenics_p_electric(my_detector) 
    my_field.fenics_p_w_electric(my_detector)
    my_field.change_data_form(my_detector)
    my_field.cal_field(my_detector)  
    ### define the tracks and type of incident particles
    my_track = Tracks()
    my_track.t_mip([50,0],[50,100],10)
    ROOT.gROOT.SetBatch(1)
    for i in range (number,number+10):
    ### drift of ionized particles
        my_drift = Drifts(my_track)
        my_drift.ionized_drift(my_track,my_field,my_detector)
        ### after the electronics
        my_electronics = Amplifier()
        qtot,ele_current=my_electronics.CSA_amp(my_detector,t_rise=0.4,t_fall=0.2,trans_imp=10)
        print("total_charge:%s fc"%(qtot*1e15))
        c = ROOT.TCanvas("Plots", "Plots", 1000, 1000)
        #ele_current.Draw("HIST")
        #c.Update()
        #c.SaveAs(output+"t_"+str(i)+"_events.C")

def main():
    args = sys.argv[1:]
    model = args[0]

    if model in ["2D","2D_SiC_PIN"]:
        twoD_time()
    
    elif model in ["2D","2D_LASER_SI_PIN"]:
        twoD_TCT()

    elif model == "2D_scan":
        output = args[1]
        v_number = int(args[2])
        twoD_time_scan(output,v_number)

    elif model == "2D_SiC_LGAD":
        twoD_time("SiC_LGAD")

    else:
        raise NameError(model) 

        
    
if __name__ == '__main__':
    #record run time
    starttime = time.time()
    main() 
    endtime=time.time()
    dtime = endtime-starttime
    print ("the process run time %.8s s" %dtime) 