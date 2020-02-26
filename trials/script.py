import subprocess as sub
import numpy as np
cf_l  = [2] #[1,2,3]
cff_l = ["land_"] #["MNIST_","land_","video_"]
cm_l  = [6] #[4,5,6,7,8]
cmm_l = ["_hosvd"] #["_mps","_cpd","_hosvd","_hooi","_random"]
#dp_l  = [["1D-2"],["2.5D-2"],["5D-2"],["1D-1"]]
for (cf,cff) in zip(cf_l,cff_l): # SELECT DATASET
    if (cf==1):
        fr_l = [[100]] #[[100],[250],[500],[1000]]
    if (cf==2):
        fr_l = [[112,240]] #[[112,240]] #[[112,240],[149,320],[298,640],[478,1024]]
    if (cf==3):
        fr_l = ["0"]
    for fr in fr_l: # SELECT RESOLUTION
        fr_str = [str(i) for i in fr]
        print("resolution: "+" ".join(fr_str))
        for (cm,cmm) in zip (cm_l,cmm_l): # SELECT COMPRESSION
            with open("../data/error_"+cff+"_".join(fr_str)+cmm+".dat", 'w') as fout:
                    fout.write("")
            dp_l  = [[np.floor(fr[0]*0.655*i).astype('i8'),np.floor(fr[1]*0.655*i).astype('i8'),3] for i in np.arange(0.1,1.1,0.1)]
            for dp in dp_l:
                dp_str = [str(i) for i in dp]
                print("decomposition parameters: "+" ".join(dp_str))
                with open('choose_file.dat', 'w') as f1:
                    f1.write(str(cf))
                with open('file_resolution.dat', 'w') as f2:
                    f2.write(" ".join(fr_str))
                with open('choose_method.dat', 'w') as f3:
                    f3.write(str(cm))
                with open('decomposition_parameters.dat', 'w') as f4:
                    f4.write(" ".join(dp_str))
                sub.call(["./main.x"])
                real_shape = fr.append(3)
                original = np.loadtxt("../data/"+cff+"_".join(fr_str)+".csv", skiprows=2)
                original = original.reshape(real_shape)
                original = original.astype(np.float64)
                original = original - np.min(original)
                original = original / np.max(original)
                original_params = np.prod(real_shape)
                reco_params = np.prod(dp)
                for real_s, rank_s in zip(fr,dp):
                    reco_params += (real_s * rank_s)
                reco = np.loadtxt("../data/"+cff+"_".join(fr_str)+cmm+"_"+"_".join(dp_str)+".csv", skiprows=0)
                reco = reco.reshape(real_shape)
                reco = reco.astype(np.float64)
                reco = reco - np.min(reco)
                reco = reco / np.max(reco)
                mse = np.mean((original-reco)**2)
                with open("../data/error_"+cff+"_".join(fr_str)+cmm+".dat", 'a') as fout:
                    fout.write(str(reco_params)+"\t")
                    fout.write(str(mse)+"\n")
