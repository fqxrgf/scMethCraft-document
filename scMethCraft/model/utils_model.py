import numpy as np
import h5py
import torch
import logging

def m_round(x):
    return int(np.round(x))


def get_logger(filename, verbosity=1, name=None):
    level_dict = {0: logging.DEBUG, 1: logging.INFO, 2: logging.WARNING}
    formatter = logging.Formatter(
        "[%(asctime)s][%(filename)s][line:%(lineno)d][%(levelname)s] %(message)s"
    )
    logger = logging.getLogger(name)
    logger.setLevel(level_dict[verbosity])

    fh = logging.FileHandler(filename, "w")
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    sh = logging.StreamHandler()
    sh.setFormatter(formatter)
    logger.addHandler(sh)

    return logger


    
def load_state(input_path,filename,load_range):

    
    if load_range != False:
        data = np.load(input_path+filename)[load_range[0]:load_range[1]]
    else:
        data = np.load(input_path+filename)
    
    return data.astype(float)



    
