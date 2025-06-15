import torch
import numpy
def output_enhanced_data(scMethCraft_part1,scMethCraft_part2,train_dataloader,cell,device):
    all_result = torch.empty(0,cell)

    scMethCraft_part1.eval()
    scMethCraft_part2.eval()

    with torch.no_grad():
        for step,[onehot, _, kmer,pos] in enumerate(train_dataloader):
            result = scMethCraft_part1(onehot.to(device).float(),kmer.to(device).float(),pos.to(device).float())
            reconstructed_matrix = torch.sigmoid(result)
            imputed_matrix = reconstructed_matrix#*(torch.isnan(targets.to(device))).int()+torch.nan_to_num(targets.to(device))
            result = scMethCraft_part2(imputed_matrix.float())
            reconstructed_matrix = torch.sigmoid(result)
            result = reconstructed_matrix#*(torch.isnan(targets.to(device))).int()+torch.nan_to_num(targets.to(device))
            all_result = torch.cat([all_result,result.cpu()])
    # all_result  = sigmoid(all_result).numpy()
    all_result  = all_result.numpy().T