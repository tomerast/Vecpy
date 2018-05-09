
import numpy as np
import scipy
from scipy.stats import norm
from scipy import io
import matplotlib.pyplot as plt
from weighted_median import *


class vec_properties:
    def __init__(self,source,ws,time_unit,time_scale_to_seconds,length_unit,length_scale_to_meter):
        self.source = source
        self.ws = ws
        self.time_unit = time_unit
        self.time_scale_to_seconds = time_scale_to_seconds
        self.length_unit = length_unit
        self.length_scale_to_meter = length_scale_to_meter
        self.velocity_units = length_unit+'/'+time_unit
        
        
    def show(self):
        print(
        'source: ',self.source,'\n',
        'window size: ',self.ws,'\n',
        'dt: ',self.dt,'\n',
        'pixel to meter: ',self.pixel_to_meter,'\n',
        'velocity x units: ',self.vxunits,'\n',
        'velocity y units: ',self.vyunits)


class field_properties:
    def __init__(self,frame,time,images_path,source,time_unit,time_scale_to_seconds,length_unit,length_scale_to_meter):
        self.frame = frame
        self.time = time
        self.images_path = images_path
        self.source = source
        self.history = ''
        self.time_unit = time_unit
        self.time_scale_to_seconds = time_scale_to_seconds
        self.length_unit = length_unit
        self.length_scale_to_meter = length_scale_to_meter
        self.velocity_units = length_unit+'/'+time_unit
    
    def show(self):
        print(
        'frame: ',self.frame,'\n',
        'absolute time: ',self.time,'\n',
        'images_path: ',self.images_path,'\n',
        'source: ',self.source,'\n',
        'dt: ',self.dt,'\n',
        'pixel to meter: ',self.pixel_to_meter,'\n',
        'x units: ',self.xunits,'\n',
        'y units: ',self.yunits,'\n',
        'velocity x units: ',self.vxunits,'\n',
        'velocity y units: ',self.vyunits)

class run_properties:
    pass
        
class vector:
    def __init__(self,X,Y,U,V,S2N,properties):
        self.X = X
        self.Y = Y
        self.U = U
        self.V = V
        self.S2N = S2N
        self.properties = properties
    
    def convert_units(self,output_length_unit,output_time_unit):
        LS = {'mm':0.001, 'cm':0.01, 'm':1.0,'meter':1.0, 'km':1000.}
        TS = {'ms':0.001, 's':1.0,'second':1.0, 'min':60.,'h':3600.,'hour':3600.}
        LS[self.properties.length_unit]=float(self.properties.length_scale_to_meter)
        TS[self.properties.time_unit]=float(self.properties.time_scale_to_seconds)
        
        self.X = self.X*(LS[self.properties.length_unit]/LS[output_length_unit])
        self.Y = self.Y*(LS[self.properties.length_unit]/LS[output_length_unit])

        self.U = self.U*(LS[self.properties.length_unit]/LS[output_length_unit])*(TS[output_time_unit]/TS[self.properties.time_unit])
        self.V = self.V*(LS[self.properties.length_unit]/LS[output_length_unit])*(TS[output_time_unit]/TS[self.properties.time_unit])
        
        self.properties.lenght_unit = output_length_unit
        self.properties.length_scale_to_meter = LS[output_length_unit]
        self.properties.time_unit = output_time_unit
        self.properties.time_scale_to_seconds = TS[output_time_unit]
        self.properties.velocity_units = output_length_unit+'/'+output_time_unit
    

class field:
    def __init__(self,field_properties):  
        self.data = {}
        self.filtered = {}
        self.properties = field_properties


    def __add__(self,other):
        check_list = []
        check_list.append(self.properties.lenght_unit == other.properties.lenght_unit)
        check_list.append(self.properties.length_scale_to_meter == other.properties.length_scale_to_meter)
        check_list.append(self.properties.time_unit == other.properties.time_unit)
        check_list.append(self.properties.time_scale_to_seconds == other.properties.time_scale_to_seconds)
        check_list.append(self.properties.velocity_units == other.properties.velocity_units)

        if all(check_list):
            sum_properties = self.properties
            sum_properties.source = 'Sum'
            sum_properties.frame = self.properties.frame + ' & ' + other.properties.frame
            sum_properties.time = self.properties.time + ' & ' + other.properties.time
            sum_properties.images_path = self.properties.images_path + ' & ' + other.properties.images_path

            sum_field = field(sum_properties)
            for xy in list(self.data.keys()):
                sum_field.add_vec(self.data[xy])
            for xy in list(other.data.keys()):
                sum_field.add_vec(other.data[xy])
    
            return sum_field

        else:
            print( 'Field properties do not match')

    def add_vec(self, vector):
        self.data[vector.X,vector.Y] = vector

    def check_if_grid_point_exists(self,x,y):
        xy = list(self.data.keys())
        return (x,y) in xy
        
    def move_to_filtered(self,vector):
        self.filtered[vector.X,vector.Y] = vector
        self.remove_vec(0,0,vector)

    def transfer(self,other):
        for xy in list(other.data.keys()):
            self.add_vec(other.data[xy])
    
    def convert_field_units(self,output_length_unit,output_time_unit):    
        XY = list(self.data.keys())
        for xy in XY:
            self.data[xy].convert_units(output_length_unit,output_time_unit)
            
        self.properties.lenght_unit = self.data[XY[0]].properties.lenght_unit
        self.properties.length_scale_to_meter = self.data[XY[0]].properties.length_scale_to_meter
        self.properties.time_unit = self.data[XY[0]].properties.time_unit
        self.properties.time_scale_to_seconds = self.data[XY[0]].properties.time_scale_to_seconds
        self.properties.velocity_units = self.data[XY[0]].properties.velocity_units
        
    def remove_vec(self,X,Y,vector=None):
        if vector is not None:
            del self.data[vector.X,vector.Y]
        else:
            del self.data[X,Y]
    
    def return_vel(self,x,y):
        u = self.data[x,y].U
        v = self.data[x,y].V

        return u,v

    def return_n_closest_neighbors(self,x,y,n=4):
        X,Y = self.return_grid()
        dist = np.sqrt((X-x)**2+(Y-y)**2)
        n_closest_neighbors = [ [(X[ind],Y[ind]),dist[ind]] for ind in dist.argsort()[:n]]
        
        return n_closest_neighbors

    def return_closest_neighbors_radius(self,x,y,radius):
        X,Y = self.return_grid()
        dist = np.sqrt((X-x)**2+(Y-y)**2)
        indecies = np.where(dist<radius)
        closest_neighbors = [[(X[indecies[0][i]],Y[indecies[0][i]]),dist[indecies[0][i]]] for i in range(len(indecies[0]))]
        return closest_neighbors

    def return_grid(self):
        XY = list(self.data.keys())
        X,Y = zip(*XY)
        X = np.array(X)
        Y = np.array(Y)
        return X,Y
    
    def return_all_velocities(self):
        XY = list(self.data.keys())
        U = np.array([self.data[xy[0],xy[1]].U for xy in XY])
        V = np.array([self.data[xy[0],xy[1]].V for xy in XY])
        return U,V

    def sub_average(self):
        XY = list(self.data.keys())
        umean,ustd,vmean,vstd = self.mean_velocity()
        for i in range(len(XY)):
            self.data[XY[i]].U = self.data[XY[i]].U - umean
            self.data[XY[i]].V = self.data[XY[i]].V - vmean

    def create_mesh_grid(self):
        X,Y = self.return_grid()
        U,V = self.return_all_velocities()
        X_mesh_grid = sorted(list(set(X)))
        Y_mesh_grid = sorted(list(set(Y)))
        X_mesh_grid,Y_mesh_grid = np.meshgrid(X_mesh_grid,Y_mesh_grid)
        U_mesh_grid = np.empty(X_mesh_grid.shape)
        U_mesh_grid.fill(np.nan)
        V_mesh_grid = np.empty(X_mesh_grid.shape)
        V_mesh_grid.fill(np.nan)
        for vec_ind in range(len(X)):
            x = X[vec_ind]
            y = Y[vec_ind]
            col = np.array(np.where(X_mesh_grid[0,:]==x))[0,0]
            row = np.array(np.where(Y_mesh_grid[:,0]==y))[0,0]
            U_mesh_grid[row,col] = U[vec_ind]
            V_mesh_grid[row,col] = V[vec_ind]
        
        return X_mesh_grid,Y_mesh_grid[::-1],U_mesh_grid[::-1],V_mesh_grid[::-1]

    def s2n_filter(self,threshold):
        XY = list(self.data.keys())
        for xy in XY:
            if self.data[xy].S2N < threshold:
                self.move_to_filtered(self.data[xy])
        
    def hist_filter(self,percentage):
        hist_u,hist_V,hist2d = self.velocity_histogram()
        XY = list(self.data.keys())
        number_of_vectors = len(XY)
        for xy in XY:
            u = self.data[xy].U
            v = self.data[xy].V
            #strech boundry edges
            hist2d[1][0] = hist2d[1][0]-1
            hist2d[1][-1] = hist2d[1][-1]+1
            hist2d[2][0] = hist2d[2][0]-1
            hist2d[2][-1] = hist2d[2][-1]+1
            
            U_histpos = np.digitize(u,hist2d[1])-1
            V_histpos = np.digitize(v,hist2d[2])-1
            if hist2d[0][U_histpos,V_histpos] / number_of_vectors < percentage/100:
                self.move_to_filtered(self.data[xy])
    
    def Z_filter(self,threshold,neighbors=4,power=1):
        XY = list(self.data.keys())
        for xy in XY:
            u = self.data[xy].U
            v = self.data[xy].V
            closest_neighbors = self.return_n_closest_neighbors(self.data[xy].X,self.data[xy].Y,neighbors+1)[1:]
            neighbor_pos , dis = zip(*closest_neighbors)
            weights = [(1/d)**power for d in dis]
            U,V = zip(*[self.return_vel(pos[0],pos[1]) for pos in neighbor_pos])

            median_U = weighted_median(U,weights)
            median_V = weighted_median(V,weights)
            median_absolute_deviation_U = weighted_median([np.abs(u_neighbor - median_U) for u_neighbor in U],weights)
            median_absolute_deviation_V = weighted_median([np.abs(v_neighbor - median_V) for v_neighbor in V],weights)

            if 0.6745*(u - median_U) / max(median_absolute_deviation_U,0.01) > threshold:
                self.move_to_filtered(self.data[xy])
                continue

            if 0.6745*(v - median_V) / max(median_absolute_deviation_V,0.01) > threshold:
                self.move_to_filtered(self.data[xy])
                continue

            
    def mean_velocity(self):
        U,V = self.return_all_velocities()
        return np.nanmean(U),np.nanstd(U),np.nanmean(V),np.nanstd(V)
    
    def velocity_histogram(self,bins=10):
        U,V = self.return_all_velocities()
        hist_U = np.histogram(U,bins)
        hist_V = np.histogram(V,bins)
        hist2d = np.histogram2d(U, V, bins)
        return hist_U,hist_V,hist2d
    
    def extract_area(self,x_boundry,y_boundry):
        area = field(self.properties)
        X,Y = self.return_grid()
        for i in range(len(X)):
            if x_boundry[0]<=X[i]<=x_boundry[1] and y_boundry[0]<=Y[i]<=y_boundry[1]:
                area.add_vec(self.data[X[i],Y[i]])

        return area

    def vel_gradients(self):
        X,Y,U,V = self.create_mesh_grid()
        Udx,Udy = np.gradient(U)
        Vdx,Vdy = np.gradient(V)
        
        return X,Y,Udx,Udy,Vdx,Vdy

    def profile(self,axis='y'):
        X,Y,U,V = self.create_mesh_grid()
        if axis=='y' or axis=='Y':
            U_profile = np.nanmean(U,axis=1)[::-1]
            V_profile = np.nanmean(V,axis=1)[::-1]
            Y_profile = Y[:,0]

            return U_profile,V_profile,Y_profile
        else: 
            U_profile = np.nanmean(U,axis=0)[::-1]
            V_profile = np.nanmean(V,axis=0)[::-1]
            X_profile = X[0,:]

            return U_profile,V_profile,X_profile
    
    def inverse_distance_interpolation(self,x,y,number_of_neighbors=4,radius=None,inverse_power=2):
        if radius is not None:
            indecies,distances = zip(*self.return_closest_neighbors_radius(x,y,radius))
        else:
            indecies,distances = zip(*self.return_n_closest_neighbors(x,y,n=4))
        weights = np.array(distances)**-float(inverse_power)
        neighbors_vel = [self.return_vel(ind[0],ind[1]) for ind in indecies]
        inter_u = (np.multiply(np.array([v[0] for v in neighbors_vel]),weights)).sum() / weights.sum()
        inter_v = (np.multiply(np.array([v[1] for v in neighbors_vel]),weights)).sum() / weights.sum()

        return inter_u,inter_v
    
    def remap(self,X,Y):
        new_feild = field(self.properties)
        X = X.flatten()
        Y = Y.flatten()
        Xold,Yold = self.return_grid()
        vec_properties = self.data[Xold[0],Yold[0]].properties
        vec_properties.source = 'Interpolation'
        for ind in range(len(X)):
            u,v = self.inverse_distance_interpolation(X[ind],Y[ind])
            vec = vector(X[ind],Y[ind],u,v,0,vec_properties)
            new_feild.add_vec(vec)
        
        self.filtered = self.data
        self.data = {}
        self.transfer(new_feild)
        
        
    
class run:
    def __init__(self,run_properties):  
        self.fields = {}
        self.properties = run_properties
    
    def add_field(self,field):
        self.fields[field.properties.frame] = field
    
    def frames(self):
        return list(self.fields.keys())
        
    def remove_field(self,frame,field=None):
        if field is not None:
            del self.fields[field.properties.frame]
        else:
            del self.fields[frame]
    
    def gp_exists_all_frames(self,x,y):
        frames = self.frames()
        gp_exists = [self.fields[f].check_if_grid_point_exists(x,y) for f in frames]
        if all(gp_exists):
            return True
        else:
            no_gp_frames = [x for x, y in zip(frames, gp_exists) if y == False]
            frames_with_gp = [x for x, y in zip(frames, gp_exists) if y == True]
            print('Frames without the requested grid point ','(',x,',',y,')',': ',no_gp_frames)
            return frames_with_gp

    def run_grid(self):
        frames = self.frames()
        Y_agp = []
        X_agp =[]
        for frame in frames:
            X,Y = self.fields[frame].return_grid()
            Y_agp += Y.tolist()
            Y_agp = sorted(list(set(Y_agp)))
            X_agp += X.tolist()
            X_agp = sorted(list(set(X_agp)))
        
        return np.meshgrid(np.array(X_agp),np.array(Y_agp))

    def grid_point_velocity(self,x,y,frames=None):
        if frames==None:
            frames = self.frames()    
            if self.gp_exists_all_frames(x,y):
                U = []
                V = []
                for f in frames:
                    u,v = self.fields[f].return_vel(x,y)
                    U.append(u)
                    V.append(v)
                    
                U = np.array(U)
                V = np.array(V)
                return U,V
        else:
            U = []
            V = []
            for f in frames:
                u,v = self.fields[f].return_vel(x,y)
                U.append(u)
                V.append(v)
                
            U = np.array(U)
            V = np.array(V)
            return U,V

    def mean_gp_velocity(self,x,y):
        for_all_frames = self.gp_exists_all_frames(x,y)
        if for_all_frames==True:
            U,V = self.grid_point_velocity(x,y)
            U_rms = U - np.nanmean(U)
            V_rms = V - np.nanmean(V)
            return np.nanmean(U),U_rms,np.nanmean(V),V_rms
        else:
            U,V = self.grid_point_velocity(x,y,for_all_frames)
            U_rms = U - np.nanmean(U)
            V_rms = V - np.nanmean(V)
            return np.nanmean(U),U_rms,np.nanmean(V),V_rms

    def mean_velocity_properties(self):
        frames = self.frames()
        U_mean = []
        V_mean = []
        for f in frames:
            u_mean,u_std,v_mean,v_std = self.fields[f].mean_velocity()
            U_mean.append(u_mean)
            V_mean.append(v_mean)
        Um = np.mean(U_mean)
        Vm = np.mean(V_mean)
        U_rms = [(np.sqrt((u-Um)**2)) for u in U_mean]
        V_rms = [(np.sqrt((v-Vm)**2)) for v in V_mean]
        print('Max in mean U velocity accures in frame: ',frames[U_mean.index(max(U_mean))])
        print('Max in mean V velocity accures in frame: ',frames[V_mean.index(max(V_mean))])
        U_mean = np.array(U_mean)
        V_mean = np.array(V_mean)
        U_rms = np.array(U_rms)
        V_rms = np.array(V_rms)
        return U_mean,U_rms,V_mean,V_rms

    def run_mean_velocities(self):
        X,Y = self.run_grid()
        U_mean = np.zeros(X.shape)
        V_mean = np.zeros(Y.shape)
        for row in range(X.shape[0]):
            for col in range(X.shape[0]):
                u,urms,v,vrms = self.mean_gp_velocity(X[row,col],Y[row,col])
                U_mean[row,col] = u
                V_mean[row,col] = v

        return U_mean,V_mean

    def mean_profile(self,axis='y'):
        frames = self.frames()
        if axis=='y' or axis=='Y':
            Y_agp = []
            for frame in frames:
                X,Y = self.fields[frame].return_grid()
                Y_agp += Y.tolist()
                Y_agp = sorted(list(set(Y_agp)))
                
            U_profiles = np.empty((len(Y_agp),len(frames)))
            U_profiles.fill(np.nan)
            V_profiles = np.empty((len(Y_agp),len(frames)))
            V_profiles.fill(np.nan)
            Y_agp = np.array(Y_agp)
            
            for col_ind,frame in list(enumerate(frames)):
                U_cur_prof,V_cur_prof,Y_cur_prof = self.fields[frame].profile(axis=axis)
                for i in range(len(Y_cur_prof)):
                     row_ind = np.where(Y_agp==Y_cur_prof[i])[0][0]
                     U_profiles[row_ind,col_ind] = U_cur_prof[i]
                     V_profiles[row_ind,col_ind] = V_cur_prof[i]
                     
            U_mean_profile = np.nanmean(U_profiles,axis=1)[::-1]
            U_number_of_vectors = np.sum(np.invert(np.isnan(U_profiles)),1)
            V_mean_profile = np.nanmean(V_profiles,axis=1)[::-1]
            V_number_of_vectors = np.sum(np.invert(np.isnan(V_profiles)),1)
            
            return U_mean_profile,U_number_of_vectors,V_mean_profile,V_number_of_vectors,Y_agp
        
        else:
            X_agp = []
            for frame in frames:
                X,Y = self.fields[frame].return_grid()
                X_agp += X.tolist()
                X_agp = sorted(list(set(X_agp)))
                
            U_profiles = np.empty((len(X_agp),len(frames)))
            U_profiles.fill(np.nan)
            V_profiles = np.empty((len(X_agp),len(frames)))
            V_profiles.fill(np.nan)
            X_agp = np.array(X_agp)
            
            for col_ind,frame in list(enumerate(frames)):
                U_cur_prof,V_cur_prof,X_cur_prof = self.fields[frame].profile(axis=axis)
                for i in range(len(X_cur_prof)):
                     row_ind = np.where(X_agp==X_cur_prof[i])[0][0]
                     U_profiles[row_ind,col_ind] = U_cur_prof[i]
                     V_profiles[row_ind,col_ind] = V_cur_prof[i]
                     
            U_mean_profile = np.nanmean(U_profiles,axis=1)[::-1]
            U_number_of_vectors = np.sum(np.invert(np.isnan(U_profiles)),1)
            V_mean_profile = np.nanmean(V_profiles,axis=1)[::-1]
            V_number_of_vectors = np.sum(np.invert(np.isnan(V_profiles)),1)
            
            return U_mean_profile,U_number_of_vectors,V_mean_profile,V_number_of_vectors,X_agp
        

    def corr_s(self,x1,y1,x2,y2):
        if self.gp_exists_all_frames(x1,y1) and self.gp_exists_all_frames(x2,y2):
            Umean1,U1,Vmean1,V1 = self.mean_gp_velocity(x1,y1)
            Umean2,U2,Vmean2,V2 = self.mean_gp_velocity(x2,y2)

            return np.inner(U1, U2)/(np.sqrt(np.inner(U1, U1)*np.inner(U2, U2))) , np.inner(V1, V2)/(np.sqrt(np.inner(V1, V1)*np.inner(V2, V2)))

        else:
            frames1 = self.gp_exists_all_frames(x1,y1)
            frames2 = self.gp_exists_all_frames(x2,y2)
            frames_for_both = set(frames1).intersection(frames2)
            U1,V1 = self.grid_point_velocity(x1,y1,frames_for_both)
            U2,V2 = self.grid_point_velocity(x2,y2,frames_for_both)
            Umean1,nouse1,Vmean1,nouse2 = self.mean_gp_velocity(x1,y1)
            Umean2,nouse1,Vmean2,nouse2 = self.mean_gp_velocity(x2,y2)
            U1 = U1 - Umean1
            V1 = V1 - Vmean1
            U2 = U2 - Umean2
            V2 = V2 - Vmean2

            return np.inner(U1, U2)/(np.sqrt(np.inner(U1, U1)*np.inner(U2, U2))) , np.inner(V1, V2)/(np.sqrt(np.inner(V1, V1)*np.inner(V2, V2)))
    
    def corr_t(self,x,y,dframes,tau=None):
        if tau is not None:
            dframes = int(tau//float(self.fields[self.frames()[0]].properties.time_scale_to_seconds))
        if self.gp_exists_all_frames(x,y):
            Umean1,U1,Vmean1,V1 = self.mean_gp_velocity(x,y)
            U2 = U1[dframes:]
            V2 = V1[dframes:]
            U1 = U1[:-dframes]
            V1 = V1[:-dframes]
            return np.inner(U1, U2)/(np.sqrt(np.inner(U1, U1)*np.inner(U2, U2))) , np.inner(V1, V2)/(np.sqrt(np.inner(V1, V1)*np.inner(V2, V2)))

        else:
            frames = self.gp_exists_all_frames(x,y)
            U1,V1 = self.grid_point_velocity(x,y,frames)
            U1 = U1 - np.nanmean(U1)
            V1 = V1 - np.nanmean(V1)
            U2 = U1[dframes:]
            V2 = V1[dframes:]
            U1 = U1[:-dframes]
            V1 = V1[:-dframes]
            return np.inner(U1, U2)/(np.sqrt(np.inner(U1, U1)*np.inner(U2, U2))) , np.inner(V1, V2)/(np.sqrt(np.inner(V1, V1)*np.inner(V2, V2)))

    def tke_frame(self,frame):
        X,Y = self.fields[frame].return_grid()
        U,V = self.fields[frame].return_all_velocities()
        tke = np.zeros(U.shape)
        for ind in range(len(X)):
            umean,urms,vmean,vrms =  self.mean_gp_velocity(X[ind],Y[ind])
            # need to check equation (devide by mean velocities???)
            tke[ind] = 0.5*(np.sqrt((U[ind]-umean)**2+(V[ind]-vmean)**2))

        return X,Y,tke

    def run_tke(self):
        X,Y = self.run_grid()
        tke = np.zeros(X.shape)
        for row in range(X.shape[0]):
            for col in range(X.shape[0]):
                u,urms,v,vrms = self.mean_gp_velocity(X[row,col],Y[row,col])
                # need to check equation (devide by mean velocities???)
                tke[row,col] = np.sqrt((1/len(urms))*np.sum(urms**2+vrms**2))

        return X,Y,tke



#debug
'''
frame = 1
dt = 0.001
length_scale = 'pixel'
length_to_meter = 0.00011946

def create_syn_prop(dt,length_scale,length_to_meter):
    global frame
    vprop = vec_properties('Tomers algorithm','64x64','dt',str(dt),str(length_scale),str(length_to_meter))
    cur_frame = (5-len(str(frame)))*str(0)+str(frame)
    fprop = field_properties(cur_frame,frame*dt,'c:/user/user/....','Tomers algorithm','dt',str(dt),str(length_scale),str(length_to_meter))
    frame+=1
    return fprop,vprop


def create_syn_field(rows,cols,xboundries,yboundries,dt,length_scale,length_to_meter,rand=None,number_of_vectors=None):
    fprop,vprop = create_syn_prop(dt,length_scale,length_to_meter)
    field1 = field(fprop)
    if rand==True:
        for i in range(number_of_vectors):
            pixel_loc = np.random.randint(0,1920,2)
            U = 5*(np.random.rand()-0.5)
            V = 1*(np.random.rand()-0.5)
            vec = vector(pixel_loc[0],pixel_loc[1],U,V,0,vprop)
            field1.add_vec(vec)
        return field1
    else:
        U = np.zeros((rows,cols))
        V = np.zeros((rows,cols))
        X,Y = np.meshgrid(np.linspace(xboundries[0],xboundries[1],cols),np.linspace(yboundries[0],yboundries[1],rows))
        for row in range(rows):
            for col in range(cols):
                u = 5*(np.random.rand()-0.5)
                v = 1*(np.random.rand()-0.5)
                vec = vector(X[row,col],Y[row,col],u,v,0,vprop)
                field1.add_vec(vec)
        return field1

def create_syn_run(number_of_fields,rows,cols,xboundries,yboundries,dt,length_scale,length_to_meter,rand=None,number_of_vectors=None):
    run1 = run('a')
    for i in range(number_of_fields):
        field1 = create_syn_field(rows,cols,xboundries,yboundries,dt,length_scale,length_to_meter,rand=None,number_of_vectors=None)
        run1.add_field(field1)
    return run1


run1 = create_syn_run(200,30,30,[0,1920],[0,1920],dt,length_scale,length_to_meter,rand=None,number_of_vectors=None)
run1.corr_t(0,0,0,0.005)        
        
'''