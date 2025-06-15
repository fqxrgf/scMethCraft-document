import torch
import pandas as pd


p_max_map_human = dict(pd.DataFrame([[1, 2,3, 4, 5, 6, 7, 8, 9, 10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25],[248956422, 242193529, 198295559, 190214555,
181538259, 170805979, 159345973, 145138636, 138394717,
133797422, 135086622, 133275309, 114364328, 107043718,
101991189, 90338345, 83257441, 80373285, 58617616, 64444167,
46709983, 50818468, 156040895, 572274151,10000]]).T.values)

p_max_map_mouse = dict(pd.DataFrame([['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',  
        '11', '12', '13', '14', '15', '16', '17', '18', '19','23', '24','25'],[195471971, 182113224, 160039680, 156508116,
              151834684, 149736546, 145441459, 129401213, 124595110,
              130694993, 122082543, 120129022, 120421639, 124902244,
              104043685, 98207768, 94987271, 90702639, 61431566, 171031299, 91744698,10000]]).astype(int).T.values)


def f(i,d_emb):
    return 1e-4+i*(d_emb/2-1-1e-4)/(d_emb/2-1)

def PE(j,d_emb,w,t):
    if j == 0:
        return t
    if j == 1:
        return t
    if j % 2 == 0:
        return torch.cos(f((j / 2 -1 ),d_emb) * w)
    if j % 2 == 1:
        return torch.sin(f(((j-1) / 2 -1 ),d_emb) * w)
    
def return_pos(pos,species = "human"):
    if species == "human":
        p_max_map = p_max_map_human
    elif species == "mouse":
        p_max_map = p_max_map_mouse
    p_max = torch.tensor([p_max_map[i] for i in pos[:,0].numpy()])
    p = (pos[:,1]+pos[:,2])/2
    t = (p/p_max).unsqueeze(1)
    w = 2*3.14159*t
    d_emb = 64
    return(torch.concatenate([PE(i,d_emb,w,t) for i in range(d_emb)],axis = 1))