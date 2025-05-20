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

def lumi_mask(tree):
    runinfo = tree["RELBO"].array()
    #print("run_number \n", runinfo[:, 0])
    #print("lumi block \n", runinfo[:, 2])
    run_number = runinfo[:, 0]
    lumi_block = runinfo[:, 2]

    keep_events = ~((run_number == 319337) & (lumi_block == 48)) ##keep good lumi blocks only
    return keep_events
  
def compute_event_mask(muon_eta, muon_phi, muon_vert, trig_eta, trig_phi):
    # Output: per-event mask (True = keep, False = reject)
    n_events = len(muon_eta)
    event_mask = []

    for evt in range(n_events):
        reject_event = False

        for i in range(len(muon_eta[evt])):
            if muon_vert[evt][i] < 0:
                continue

            for j in range(len(trig_eta[evt])):
                dR = deltaR(muon_eta[evt][i], muon_phi[evt][i],
                            trig_eta[evt][j], trig_phi[evt][j])
                if dR < 0.01:
                    reject_event = True
                    break
            if reject_event:
                break

        event_mask.append(not reject_event)  # True = keep

    return ak.Array(event_mask)

    
#def muon_trigger_object_reject(sumpt, trig_eta, trig_phi, muon_eta, muon_phi, muon_vert):
#    final_vertex_mask = []
#
#    for evt_idx in range(len(muon_eta)):
##        ##some muon vertex index were empty
#        if len(muon_vert[evt_idx]) == 0 or ak.max(muon_vert[evt_idx]) < 0:
#            final_vertex_mask.append([])
#            continue
#
#        n_verts = int(ak.max(muon_vert[evt_idx]) + 1)
#        evt_mask = [False] * n_verts
#
#        for mu_eta, mu_phi, mu_vtx in zip(muon_eta[evt_idx], muon_phi[evt_idx], muon_vert[evt_idx]):
#            if mu_vtx < 0 or mu_vtx >= n_verts:
#                continue  # skip invalid or out-of-bounds indices
#
#            for trig_eta_val, trig_phi_val in zip(trig_eta[evt_idx], trig_phi[evt_idx]):
#                dR = deltaR(mu_eta, mu_phi, trig_eta_val, trig_phi_val)
#                if dR < 0.01:
#                    evt_mask[mu_vtx] = True
#                    break  # No need to check other triggers for this muon
#            final_vertex_mask.append(evt_mask)
#
#    print(len(final_vertex_mask[0]), len(muon_vert[0]))
#
#    return ak.Array(final_vertex_mask)


#def muon_trigger_object_reject(sumpt, trig_eta, trig_phi, muon_eta, muon_phi, muon_vert):
#
#    ##need to use this but could not find the awkward version
#    #muons = vector.awkward.zip({
#    #"pt": ak.ones_like(muon_eta),
#    #"eta": muon_eta,
#    #"phi": muon_phi,
#    #"mass": ak.ones_like(muon_eta)
#    #},with_name="LorentzVector")
#
#    #trigs = vector.awkward.zip({
#    #"pt": ak.ones_like(trig_eta),
#    #"eta": trig_eta,
#    #"phi": trig_phi,
#    #"mass": ak.ones_like(trig_eta)
#    #}, with_name="LorentzVector")
#
#
#    #cartesian_muon_trig = ak.cartesian([muon_eta, muon_phi, trig_eta, trig_phi], nested=True)
#    final_vertex_mask = ak.zeros_like(sumpt, dtype=bool)
#
#    print(final_vertex_mask, len(final_vertex_mask))
#    for i in range(len(muon_eta)):
#        for j in range(len(trig_eta)):
#
#            deltaR_values = deltaR(muon_eta[i], muon_phi[i], trig_eta[j], trig_phi[j])
#            print(deltaR_values[i])
#        
#            if deltaR_values[i] < 0.01:
#                final_vertex_mask[0][i] = True
#                print(final_vertex_mask)
#
#                break
#
#    print(final_vertex_mask)
#    #muon_trig_reject = ak.any(deltaR_values < 0.01)
#
#    return vertex_reject_mask

def compute_muon_trig_reject(tree):

    muon_pt = tree["VLMuonPt"].array()
    muon_eta = tree["VLMuonEta"].array()
    muon_phi = tree["VLMuonPhi"].array()
    muon_vert = tree["MuonVertInd"].array()
    sumpt = tree["SumPt"].array()

    #trig_collections = [
    #    "Mu8p5IP3p5Part012345", "Mu10p5IP3p5Part012345", "Mu9IP6Part012345", "Mu8IP3Part012345", "Mu12IP6Part01234",
    #    "Mu9IP5Part01234", "Mu7IP4Part01234", "Mu9IP4Part01234", "Mu8IP5Part01234", "Mu8IP6Part01234", "Mu9IP0Part0", "Mu9IP3Part0"
    #]

    ##create an n-dimensional array based on decision from each trigger. then ak.sum the decisions and save it per vertex and then reject that vertex.

    #trig_collections = ["Mu8IP3Part012345"]

    event_keep_mask = ak.Array([True] * len(muon_eta))  # Start by keeping all events

    trig_collections = [
        "Mu8p5IP3p5Part012345", "Mu10p5IP3p5Part012345", "Mu9IP6Part012345", "Mu8IP3Part012345", "Mu12IP6Part01234",
        "Mu9IP5Part01234", "Mu7IP4Part01234", "Mu9IP4Part01234", "Mu8IP5Part01234", "Mu8IP6Part01234", "Mu9IP0Part0", "Mu9IP3Part0"
    ]

    for name in trig_collections:
        trig_eta = tree[f"HLT{name}Eta"].array()
        trig_phi = tree[f"HLT{name}Phi"].array()

        if(ak.sum(ak.num(trig_eta)) == 0) :
            print("No trigger objects found ")
            continue
        updated_mask = compute_event_mask(muon_eta, muon_phi, muon_vert, trig_eta, trig_phi)

        event_keep_mask = event_keep_mask & updated_mask
        #print((event_keep_mask))


    #print("event_keep_mask", (event_keep_mask))    
    return event_keep_mask


def process_file(filename, outdir="/cms/kaur/output/lumi_mask/"):
    tree = uproot.open(filename)["rootTupleTreeVeryLoose/tree"]

    lumi_select = lumi_mask(tree)
    print(lumi_select)
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

    # apply *all* track cuts
    #track_cut_names = [
    #        "Track_SumPt_min",
    #        "Track_VSP_ratio_min"
    #        ]

    #track_masks = [apply_cut(tree, name) for name in track_cut_names]
    #combined_track_masks = track_masks[0]

    #for mask in track_masks[1:]:
    #    combined_track_masks = combined_track_masks & mask  # logical AND per-vertex
    trig_select = compute_muon_trig_reject(tree)

    #print("accepted_vertices ", accepted_vertices)
    final_selection = accepted_vertices & trig_select & lumi_select 
    #print("final_selection ", final_selection)
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
