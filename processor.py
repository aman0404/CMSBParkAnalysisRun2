import yaml, uproot 
import awkward as ak
from ROOT import TH1F
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import os
import vector 
vector.register_awkward()


# 1) load all cut definitions
cfg = yaml.safe_load(open("selections.yaml"))["cuts"]


# Function to calculate delta R
def deltaR(eta1, phi1, eta2, phi2):
    delta_eta = eta1 - eta2
    delta_phi = phi1 - phi2
    delta_phi = np.arctan2(np.sin(delta_phi), np.cos(delta_phi))  # Ensure periodicity of phi
    return np.sqrt(delta_eta**2 + delta_phi**2)


def apply_cut(tree, cut_name):
    cut = cfg[cut_name]
    # fetch arrays
    if "branch2" in cut:
        a1 = tree[cut["branch1"]].array()
        a2 = tree[cut["branch2"]].array()
        rhs = eval(cut.get("expr", str(cut["value"])), {"branch2": a2})
        mask = eval(f"a1 {cut['op']} rhs")
    else:
        a = tree[cut["branch"]].array()
        mask = eval(f"a {cut['op']} {cut['value']}")
    return mask 
    #return mask if cut.get("reject_if", False) else mask


def muon_trigger_object_reject(trig_pt, trig_eta, trig_phi, muon_pt, muon_eta, muon_phi, muon_vert):

    ##need to use this but could not find the awkward version
    #muons = vector.awkward.zip({
    #"pt": ak.ones_like(muon_eta),
    #"eta": muon_eta,
    #"phi": muon_phi,
    #"mass": ak.ones_like(muon_eta)
    #},with_name="LorentzVector")

    #trigs = vector.awkward.zip({
    #"pt": ak.ones_like(trig_eta),
    #"eta": trig_eta,
    #"phi": trig_phi,
    #"mass": ak.ones_like(trig_eta)
    #}, with_name="LorentzVector")


    cartesian_muon_trig = ak.cartesian([muon_eta, muon_phi, trig_eta, trig_phi], nested=True)


    # Extract the arrays for deltaR calculation
    muon_eta_flat = cartesian_muon_trig["0"]  # First muon eta
    muon_phi_flat = cartesian_muon_trig["1"]  # First muon phi
    trig_eta_flat = cartesian_muon_trig["2"]  # Trigger eta
    trig_phi_flat = cartesian_muon_trig["3"]  # Trigger phi

    deltaR_values = deltaR(muon_eta_flat, muon_phi_flat, trig_eta_flat, trig_phi_flat)

    muon_trig_reject = deltaR_values < 0.01
    #muon_trig_reject = ak.any(deltaR_values < 0.01)

    return muon_trig_reject

def compute_muon_trig_reject(tree):

    muon_pt = tree["VLMuonPt"].array()
    muon_eta = tree["VLMuonEta"].array()
    muon_phi = tree["VLMuonPhi"].array()
    muon_vert = tree["MuonVertInd"].array()

    num_vertices = len(tree["SumPt"].array())  # assumes SumPt is per-vertex

    muon_trig_reject = ak.Array([False] * num_vertices)

    #trig_collections = [
    #    "Mu8p5IP3p5Part012345", "Mu10p5IP3p5Part012345", "Mu9IP6Part012345", "Mu8IP3Part012345", "Mu12IP6Part01234",
    #    "Mu9IP5Part01234", "Mu7IP4Part01234", "Mu9IP4Part01234", "Mu8IP5Part01234", "Mu8IP6Part01234", "Mu9IP0Part0", "Mu9IP3Part0"
    #]

    trig_collections = ["Mu8IP3Part012345"]

    for name in trig_collections:
        trig_pt = tree[f"HLT{name}Pt"].array()
        trig_eta = tree[f"HLT{name}Eta"].array()
        trig_phi = tree[f"HLT{name}Phi"].array()

        if(ak.sum(ak.num(trig_pt)) == 0 and ak.sum(ak.num(trig_eta)) == 0 and ak.sum(ak.num(trig_phi)) == 0):
            print("No trigger objects found ")
            print(trig_pt)
            continue

        muon_trig_reject = muon_trig_reject | muon_trigger_object_reject(
            trig_pt, trig_eta, trig_phi, muon_pt, muon_eta, muon_phi, muon_vert
        )
    return muon_trig_reject


def process_file(filename, outdir="/cms/kaur/output"):
    tree = uproot.open(filename)["rootTupleTreeVeryLoose/tree"]

    # apply *all* vertex cuts
    vertex_cut_names = [
    "IsGoodRecoVertex",
    "VectorSumPt_vs_SumPt",
    "VLTrackNByVert_max",
    ]


    vertex_masks = [apply_cut(tree, name) for name in vertex_cut_names]
    combined_mask = vertex_masks[0]
    for mask in vertex_masks[1:]:
        combined_mask = combined_mask & mask  # logical AND per-vertex
    
    accepted_vertices = ak.any(combined_mask, axis=1)  # accept event if any vertex passed all cuts

    muon_trig_reject_mask = (compute_muon_trig_reject(tree))
    

    # Step-by-step reduce from [ev][vert][mu][trig][bool]
    mask_lvl1 = ak.any(muon_trig_reject_mask, axis=-1)  # remove bool level
    mask_lvl2 = ak.any(mask_lvl1, axis=-1)              # over triggers
    mask_lvl3 = ak.any(mask_lvl2, axis=-1)              # over muons
    mask_lvl4 = ak.any(mask_lvl3, axis=-1)              # over vertices
    
    # Now shape is [n_events]
    muon_trig_reject_event_mask = mask_lvl4


    # apply *all* track cuts
    #track_cut_names = [
    #        "Track_SumPt_min",
    #        "Track_VSP_ratio_min"
    #        ]

    #track_masks = [apply_cut(tree, name) for name in track_cut_names]
    #combined_track_masks = track_masks[0]

    #for mask in track_masks[1:]:
    #    combined_track_masks = combined_track_masks & mask  # logical AND per-vertex

    final_selection = accepted_vertices & muon_trig_reject_event_mask 

    # now select entries
    select_branches = ["Sphericity", "Cluster1Chi2", "Cluster2Chi2", "TotalChi2", "PVX", "PVY","PVZ", "SumPt", "SumPtSquare", "Thrust", "TransverseSphericity",
                 "TransverseThrust", "DisplacedTrackN", "IsGoodRecoVertex"]
    data = tree.arrays(filter_name=select_branches, entry_stop=None)
    # Mask and flatten each branch
    out_dict = {}
    for branch in select_branches:
        out_dict[branch] = ak.flatten(data[branch][final_selection])

    #sel_sph   = ak.flatten(data["Sphericity"][accepted_vertices])
    #sel_disTracks  = ak.flatten(data["DisplacedTrackN"][accepted_vertices])



# ...fill your histograms here using sel_sumpt, sel_thrust, etc.
#    out_dict = {
#            "Sphericity": sel_sph,
#            "DisplacedTrackN": sel_disTracks,
#            }
    # Convert to flat Pandas DataFrame
    df = pd.DataFrame({k: ak.to_numpy(v) for k, v in out_dict.items()})
#    df = pd.DataFrame({
#        "Sphericity": ak.to_numpy(sel_sph),
#        "DisplacedTrackN": ak.to_numpy(sel_disTracks),
#    })

    
    os.makedirs(outdir, exist_ok=True)

    # Save to named output file
    out_file = os.path.join(outdir, os.path.basename(filename).replace(".root", ".parquet"))
    
    # Save as Parquet
    df.to_parquet(out_file, engine="pyarrow", index=False)
    print(f"Saved: {out_file}")


if __name__=="__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("files", nargs="+")
    args = p.parse_args()
    for fn in args.files:
        print(f"Processing {fn}...")
        process_file(fn)
        print(f"Finished {fn}")
