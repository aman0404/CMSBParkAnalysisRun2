import uproot
from ROOT import TH1F, TCanvas, gRandom, TH2F
import matplotlib.pyplot as plt
import awkward as ak

path = "/cms/insta/darshan/crabdirs/Run2UL/TAG106X01/Data/2018/ParkingBPH5/Run2018A-UL2018_MiniAODv2-v1__MINIAOD/241217_193735/0000/"
file = uproot.open(path+"analysisTree_147.root")

mytree = file["rootTupleTreeVeryLoose/tree"]

### reject common vertex
good_vertex = mytree["IsGoodRecoVertex"].arrays()["IsGoodRecoVertex"] == 0

### vector sum pT
vectorSumPt = (mytree["VectorSumPt"].arrays()["VectorSumPt"]) > (0.05*mytree["SumPt"].arrays()["SumPt"]+15)

#SumPt
sumPt = mytree["SumPt"].arrays()["SumPt"] >= 25.

trackNByVert = mytree["VLTrackNByVert"].arrays()["VLTrackNByVert"] < 3

rejected_vertices = (good_vertex | vectorSumPt | sumPt | trackNByVert
    # ptbin_cut |  # Include if using PT bins
    #revertex_filter
)

# Inverse to get accepted vertices (optional)
accepted_vertices = ~rejected_vertices

### common track rejection
track_reject_sumPt = (len(mytree["SumPt"].arrays()["SumPt"]) > 0 ) & (mytree["SumPt"].arrays()["SumPt"] > 20) & ((mytree["VectorSumPt"].arrays()["VectorSumPt"])/(mytree["SumPt"].arrays()["SumPt"]) >= 0.4)

commonTrackReject = track_reject_sumPt | trackNByVert

accepted_tracks = ~commonTrackReject

load_branches = ["Cluster1Chi2", "Cluster2Chi2", "TotalChi2", "PVX", "PVY","PVZ", "SumPt", "SumPtSquare", "Thrust", "TransverseSphericity",
                 "TransverseThrust", "DisplacedTrackN", "GoodVertexN", "VertexN"]

selected_data = mytree.arrays(filter_name=load_branches)

print(selected_data)
#for names in load_branches:
#    if (names == "PVX") or (names == "PVY") or (names == "PVZ") or (names == "Thrust") or (names == "TransverseThrust") or (names == "TransverseSphericity"):
#        print(names)
#        plt.hist(ak.flatten(selected_data[names]), bins=100, range=(-1, 1))
#        plt.yscale('log')
#    else:
#        print(names)
#        plt.hist(ak.flatten(selected_data[names]), bins=100, range=(0, 100))
#    plt.xlabel(names)
#    plt.ylabel("Events")
#    plt.show()


