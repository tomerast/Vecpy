import numpy as np
import scipy
from data_structures import  *
import os
import yaml


def field_data_properties(source,ws,frame,time,images_path,time_unit,time_scale_to_seconds,length_unit,length_scale_to_meter):
    vec_prop = vec_properties(str(source),str(ws),str(time_unit),str(time_scale_to_seconds),str(length_unit),str(length_scale_to_meter))
    field_prop = field_properties(str(frame),str(time),str(images_path),str(source),str(time_unit),str(time_scale_to_seconds),str(length_unit),str(length_scale_to_meter))

    return vec_prop,field_prop

def load_vec1(x,y,u,v,s2n,vec_prop):
    vec = vector(x,y,u,v,s2n,vec_prop)
    return vec

def load_field1(X,Y,U,V,S2N,field_prop,vec_prop,source=None,ws=None,frame=None,time=None,images_path=None,time_unit=None,time_scale_to_seconds=None,length_unit=None,length_scale_to_meter=None):
    #data proccessing
    if type(X) == np.ndarray and type(Y) == np.ndarray and type(U) == np.ndarray and type(V) == np.ndarray:
        if len(X.shape)==1 and len(Y.shape)==1 and len(U.shape)==1 and len(V.shape)==1:
            pass
        elif len(X.shape)==1 and len(Y.shape)==1 and len(U.shape)==2 and len(V.shape)==2:
            U = U.flatten()
            V = V.flatten()
            S2N = S2N.flatten()
        else:
            X = X.flatten()
            Y = Y.flatten()
            U = U.flatten()
            V = V.flatten()
            S2N = S2N.flatten()

    elif type(X) == list and type(Y) == list and type(U) == np.ndarray and type(V) == np.ndarray:
        X = np.array(X)
        Y = np.array(Y)

    if field_prop is None or vec_prop is None:
        vec_prop,field_prop = field_data_properties(source,ws,frame,time,images_path,time_unit,time_scale_to_seconds,length_unit,length_scale_to_meter)


    field1 = field(field_prop)
    for ind in range(len(X)):
        vec = load_vec1(X[ind],Y[ind],U[ind],V[ind],S2N[ind],vec_prop)
        field1.add_vec(vec)
    return field1

def parse_jobfile(jobfile):
    source = 'Tomers PIV algorithm'
    ws = jobfile.piv_paramters['window_size']
    if jobfile.pair_or_set.lower() == 'pair':
        if jobfile.image_files == 1:
            frame = jobfile.frame_name
        else:
            frame = [jobfile.frame_a,jobfile.frame_b]
    elif jobfile.pair_or_set.lower() == 'set':
        frame = jobfile.frames
    
    images_path = jobfile.images_path
    time_unit = 'dt'
    length_unit = 'pixel'
    time = ''
    time_scale_to_seconds = ''
    length_scale_to_meter=''

    return field_data_properties(source,ws,frame,time,images_path,time_unit,time_scale_to_seconds,length_unit,length_scale_to_meter)
    

def check_for_yaml(file_name,path):
    if os.path.isfile(path+str(file_name).split('.')[0]+'.yaml'):
        with open(path+str(file_name).split('.')[0]+'.yaml', 'r') as stream:
            try:
                job1 = yaml.load(stream)
                return parse_jobfile(job1)
            except yaml.YAMLError as exc:
                print(exc)
    else:
        return None
    

def load_field_from_npz(file_name,path,field_prop,vec_prop):
    file_name = str(file_name)
    if file_name.endswith('.npz'):
        data  = np.load(path+file_name)
    else:
        data  = np.load(path+file_name+'.npz')
    yaml_data = check_for_yaml(file_name,path)
    
    X = data['X']
    Y = data['Y']
    U = data['U']
    V = data['V']
    S2N = data['S2N']

    if yaml_data is not None:
        field1 = load_field1(X,Y,U,V,S2N,yaml_data[1],yaml_data[0])
    else:
        field1 = load_field1(X,Y,U,V,S2N,field_prop,vec_prop)

    return field1

def load_run_from_directory(directory_path):
    run1 = run('a')
    files = [ fname for fname in os.listdir(directory_path) if fname.endswith('.npz')]
    for file in files:
        field1 = load_field_from_npz(file,directory_path,'','')
        run1.add_field(field1)
    return run1
        

    