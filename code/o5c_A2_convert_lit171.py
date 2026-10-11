"""src-literature-171 (OSF rx2zd, CC-BY-4.0; Johnson EL, Lin JJ, King-Stephens D, Weber PB, Laxer KD, Saez I, Girgis F, D'Esposito M,
Knight RT, Badre D (2023) "A rapid theta network mechanism for flexible information encoding", Nat Commun 14:2872,
doi:10.1038/s41467-023-38574-7) -> iEEG-BIDS DERIVATIVE (epoched, bipolar-referenced, preprocessed iEEG).
Source: files/data/S<N>/S<N>.mat (MAT v7.3 FieldTrip struct `data`: fsample, label (bipolar pairs), elec.chanpos/elecpos [3 x n],
time{}, trial{} (samples x ch), trialinfo (18 x ntrials cell; field names in files/data/trialinfo.mat), sampleinfo),
files/data/S<N>/yeo.mat (Yeo-7 network per channel: yeo_key, yeo_num), behavior.mat (MATLAB table; sourcedata only), scripts/.
Pattern: o5c_A1_convert_lit169.py. Usage: python o5c_A2_convert_lit171.py <open500 record dir> <bids root>"""
import glob, json, os, re, sys
import numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from o5c_A2_common import write_epochs, events_json, h5str, h5cellstr, wtsv, wjson, finalize_common, copy_code, copy_sourcedata

REC, OUT = sys.argv[1], sys.argv[2]
SRC = os.path.join(REC, "files")
report = {"recordings": []}


def tval(f, o):
    if isinstance(o, h5py.Group):
        return "struct"
    mc = o.attrs.get("MATLAB_class", b"")
    if "MATLAB_empty" in o.attrs:
        return "n/a"
    if mc == b"char":
        return h5str(f, o)
    a = np.asarray(o[()])
    if a.dtype == object or mc == b"cell":
        refs = a.ravel()
        return tval(f, f[refs[0]]) if refs.size == 1 else ";".join(str(tval(f, f[r])) for r in refs)
    if a.size == 1:
        v = float(a.ravel()[0])
        if not np.isfinite(v):
            return "n/a" if np.isnan(v) else ("Inf" if v > 0 else "-Inf")
        return int(v) if v == int(v) else round(v, 6)
    return "n/a" if a.size == 0 else ";".join(str(x) for x in a.ravel().tolist()[:20])


with h5py.File(os.path.join(SRC, "data", "trialinfo.mat"), "r") as f:
    TI_NAMES = [re.sub(r"[^A-Za-z0-9_]", "_", s) for s in h5cellstr(f, f["trialinfo"])]
TI_NAMES = [n + "_trialinfo" if n in ("onset", "duration", "sample", "trial_index", "trial_time_start") else n for n in TI_NAMES]
print("trialinfo fields:", TI_NAMES)


def pairtype(lab):
    """bipolar pair 'ROF4-ROF5' -> adjacent numbers = depth track (sEEG, deep-surface); non-adjacent = grid row (ECoG) (article's montage rule)."""
    m = re.match(r"^([A-Za-z]+)(\d+)-([A-Za-z]+)(\d+)$", lab)
    if not m or m.group(1) != m.group(3):
        return "OTHER"
    return "SEEG" if abs(int(m.group(2)) - int(m.group(4))) == 1 else "ECOG"


scans = {}
subs_dirs = sorted(glob.glob(os.path.join(SRC, "data", "S*")), key=lambda p: int(os.path.basename(p)[1:]))
for sd in subs_dirs:
    sid = int(os.path.basename(sd)[1:]); sub = f"sub-{sid:02d}"; D = os.path.join(OUT, sub, "ieeg"); stem = f"{sub}_task-wmgating"
    with h5py.File(os.path.join(sd, f"S{sid}.mat"), "r") as f:
        d = f["data"]
        fs = float(np.asarray(d["fsample"][()]).ravel()[0])
        labels = h5cellstr(f, d["label"])
        chanpos = np.asarray(d["elec"]["chanpos"][()]).reshape(3, -1).T if "elec" in d else None
        elabels = h5cellstr(f, d["elec"]["label"]) if "elec" in d else labels
        T = np.asarray(d["trial"][()]).ravel(); TM = np.asarray(d["time"][()]).ravel()
        TI = np.asarray(d["trialinfo"][()])  # (18, ntrials) refs
        assert TI.shape[0] == len(TI_NAMES), TI.shape
        ntr = T.size
        assert TM.size == ntr == TI.shape[1], (T.size, TM.size, TI.shape)
        trials, times, rows = [], [], []
        for k in range(ntr):
            x = np.asarray(f[T[k]][()])
            if x.shape[1] != len(labels) and x.shape[0] == len(labels):
                x = x.T
            t = np.asarray(f[TM[k]][()]).ravel()
            assert x.shape == (t.size, len(labels)), (x.shape, t.shape, len(labels))
            rows.append([tval(f, f[TI[j, k]]) for j in range(len(TI_NAMES))])
            trials.append(x.astype(np.float64)); times.append(t)
    yeo_key = yeo_num = None
    yp = os.path.join(sd, "yeo.mat")
    if os.path.exists(yp):
        with h5py.File(yp, "r") as f:
            yeo_key = h5cellstr(f, f["yeo_key"]); yeo_num = np.asarray(f["yeo_num"][()]).ravel().astype(int).tolist()
        if len(yeo_key) != len(labels):  # yeo.mat lists the task-responsive electrodes only, without channel names: not mappable
            report.setdefault("yeo_unmapped", {})[sub] = {"n_yeo": len(yeo_key), "n_channels": len(labels), "networks": dict(zip(*[list(map(str, u)) for u in np.unique(yeo_key, return_counts=True)]))}
            yeo_key = yeo_num = None
    types = [pairtype(l) for l in labels]
    pos = {}
    if chanpos is not None:
        assert chanpos.shape[0] == len(elabels), (chanpos.shape, len(elabels))
        pos = {l: chanpos[i] for i, l in enumerate(elabels)}
    ch_extra = []
    for i, l in enumerate(labels):
        g = re.match(r"^([A-Za-z]+)", l)
        ch_extra.append({"group": g.group(1) if g else "n/a", "reference": "bipolar (" + l + ")", "low_cutoff": 0.1, "high_cutoff": 300 if fs >= 1000 else "n/a",
                         "description": f"bipolar pair of neighbouring contacts (article montage); Yeo-7 network {yeo_key[i]} ({yeo_num[i]})" if yeo_key else "bipolar pair of neighbouring contacts (article montage)"})
    r = write_epochs(D, stem, trials, fs, labels, ["µV"] * len(labels), types, {
        "TaskName": "wmgating",
        "TaskDescription": "Single-trial working-memory gating task: three sequential stimuli (a digit as context, a letter, a symbol; context first = CF, context last = CL), then a two-alternative forced choice on the target item (article). Trials are epoched from -1 s before the first stimulus to +1 s after the response; trialinfo fields as released in events.tsv.",
        "Instructions": "n/a (see article)", "PowerLineFrequency": 60,
        "SoftwareFilters": {"AuthorsPreprocessing": {"Description": "0.1 Hz high-pass and (1 kHz data) 300 Hz low-pass FIR filters, 60 Hz harmonics removed (DFT), demeaned, epoched, artifact/epileptiform electrodes and epochs removed, bipolar montage (article). UCI data resampled from 5 kHz to 1 kHz. Values written unchanged."}},
        "HardwareFilters": "n/a", "Manufacturer": "Nihon Kohden",
        "iEEGReference": "bipolar montage between neighbouring contacts of the same anatomical structure (channel names give the pair)",
        "iEEGPlacementScheme": "Clinical ECoG grids/strips and/or sEEG depth tracks (8-, 10-, 12-channel, 5 mm spacing) for epilepsy monitoring (UCI, UCD, CPMC; article). Channel type is inferred from the pair numbering (adjacent contacts = depth track, SEEG; non-adjacent = grid row, ECOG), see README.",
        "ElectricalStimulation": False,
    }, TI_NAMES, rows, ch_extra=ch_extra, comment=f"o5c_A2_convert_lit171.py: S{sid}.mat data.trial back-to-back, float32", trial_times=times)
    events_json(os.path.join(D, stem + "_events.json"), {k: {"Description": f"data.trialinfo row '{k}' as released (names from trialinfo.mat)"} for k in TI_NAMES})
    # electrodes: released chanpos (space not stated: native or MNI, see README) + Yeo network
    erows = []
    for i, l in enumerate(labels):
        v = pos.get(l)
        erows.append([l, *(["%.4f" % c for c in v] if v is not None and np.all(np.isfinite(v)) else ["n/a"] * 3), "n/a", ch_extra[i]["group"], yeo_key[i] if yeo_key else "n/a", yeo_num[i] if yeo_num else "n/a"])
    wtsv(os.path.join(D, f"{sub}_space-Other_electrodes.tsv"), ["name", "x", "y", "z", "size", "group", "yeo7_network", "yeo7_number"], erows)
    wjson(os.path.join(D, f"{sub}_space-Other_electrodes.json"), {"SpatialReference": "n/a (release does not state whether chanpos is in the subject's native space or MNI space)", "yeo7_network": {"Description": "Yeo-7 atlas network assigned by the authors to the bipolar channel from its MNI coordinates (release yeo.mat yeo_key): VSN visual, SMN somatomotor, DAN dorsal attention, VAN ventral attention, LBN limbic, FP frontoparietal, DMN default mode"},
                                                                  "yeo7_number": {"Description": "Yeo-7 network number as released (yeo_num)"}})
    wjson(os.path.join(D, f"{sub}_space-Other_coordsystem.json"), {"iEEGCoordinateSystem": "Other", "iEEGCoordinateUnits": "mm",
          "iEEGCoordinateSystemDescription": "Coordinates of the bipolar channels as released in data.elec.chanpos (FieldTrip elec structure; chanpos == elecpos). The article localizes electrodes in native space (CT to MR co-registration in FieldTrip) and transforms them to MNI space for group display and Yeo-7 assignment; the release does not state which of the two spaces chanpos is in, so the system is recorded as Other.",
          "iEEGCoordinateProcessingDescription": "Copied from the release; no processing."})
    r.update(sub=sub, sid=sid, fs=fs, types={t: types.count(t) for t in set(types)}, yeo=dict(zip(*np.unique(yeo_key, return_counts=True))) if yeo_key else None, n_coords=sum(1 for l in labels if l in pos))
    report["recordings"].append(r)
    scans.setdefault(sub, []).append([f"ieeg/{stem}_ieeg.vhdr", "n/a"])
    print(stem, json.dumps(r, default=str)[:400], flush=True)
for sub, rws in scans.items():
    wtsv(os.path.join(OUT, sub, f"{sub}_scans.tsv"), ["filename", "acq_time"], sorted(rws))
subs = sorted(scans)
wtsv(os.path.join(OUT, "participants.tsv"), ["participant_id", "age", "sex", "source_id", "sampling_frequency", "n_trials", "n_channels"],
     [[r["sub"], "n/a", "n/a", f"S{r['sid']}", r["fs"], r["n_trials"], r["n_channels"]] for r in report["recordings"]])
wjson(os.path.join(OUT, "participants.json"), {"age": {"Description": "not in the release (article Table S1)"}, "sex": {"Description": "not in the release"}, "source_id": {"Description": "subject folder in the release"},
                                               "sampling_frequency": {"Description": "Sampling rate of the released epochs (1 kHz UCI/CPMC, 512 Hz UCD)", "Units": "Hz"}, "n_trials": {"Description": "trials in the release"}, "n_channels": {"Description": "bipolar channels in the release"}})
sd = [(p, os.path.relpath(p, SRC)) for p in sorted(glob.glob(os.path.join(SRC, "**", "*"), recursive=True)) if os.path.isfile(p)]
copy_sourcedata(sd, os.path.join(OUT, "sourcedata", "osf-rx2zd"), scrub=True)
copy_code(OUT, [__file__])
Rr = report["recordings"]; hrs = sum(r["duration_s"] for r in Rr) / 3600
finalize_common(OUT,
    name="s500-literature-171: Working-memory gating task iEEG, epoched bipolar ECoG/sEEG from 11 patients with Yeo-7 network labels (Johnson et al., 2023) (derivative)",
    dtype="derivative", license_="CC-BY-4.0",
    authors=["Elizabeth L. Johnson", "Jack J. Lin", "David King-Stephens", "Peter B. Weber", "Kenneth D. Laxer", "Ignacio Saez", "Fady Girgis", "Mark D'Esposito", "Robert T. Knight", "David Badre"],
    refs=["https://doi.org/10.1038/s41467-023-38574-7", "https://doi.org/10.17605/OSF.IO/RX2ZD", "https://osf.io/rx2zd/"],
    howto="Cite Johnson EL, Lin JJ, King-Stephens D, Weber PB, Laxer KD, Saez I, Girgis F, D'Esposito M, Knight RT, Badre D (2023) A rapid theta network mechanism for flexible information encoding. Nature Communications 14:2872, doi:10.1038/s41467-023-38574-7; and the OSF dataset https://osf.io/rx2zd/.",
    keywords=["intracranial EEG", "working memory", "input gating", "theta", "high-frequency broadband", "prefrontal cortex", "Yeo networks", "epilepsy"],
    source=[{"DOI": "doi:10.17605/OSF.IO/RX2ZD", "URL": "https://osf.io/rx2zd/", "Version": "OSF node rx2zd, 44 files as listed 2026-10-11 (sha256 per file in sourcedata/osf-rx2zd/DEIDENTIFICATION_MANIFEST.tsv)"}],
    generated_by=[{"Name": "o5c_A2_convert_lit171.py", "Description": "Repackaging of the FieldTrip epoch structs (data.trial) into BrainVision float32, trials back-to-back, with data.trialinfo as events, data.elec.chanpos as electrodes and the Yeo-7 assignment as electrode/channel labels. Code in code/."}],
    ethics=["Collected under research protocols approved by the Institutional Review Boards of UC Irvine, UC Davis and California Pacific Medical Center; written informed consent (article)."],
    readme=f"""# Johnson et al. (2023), "A rapid theta network mechanism for flexible information encoding": epoched iEEG

Preprocessed, bipolar-referenced intracranial EEG from {len(subs)} neurosurgical patients ({len(Rr)} recordings, {hrs:.2f} h of
epochs) performing a single-trial working-memory gating task (context-first / context-last), released by the authors on OSF
(https://osf.io/rx2zd/, CC-BY 4.0) together with the analysis scripts.

THIS IS A DERIVATIVE DATASET: the release holds the authors' epoched, filtered, artifact-cleaned bipolar channels (article
Methods), not the raw clinical recordings.

## Contents
- `sub-NN/ieeg/sub-NN_task-wmgating_ieeg.*`: the released `data.trial` epochs written back-to-back (one `New Segment` marker
  per trial; `events.tsv` gives each trial's position, its authors' time axis start (-1 s) and every `trialinfo` field as
  released: {', '.join(TI_NAMES)}). Unit µV assumed (not stated in the release); values unchanged.
- Channel names are the bipolar pairs (e.g. `LF9-LF1`). The channel type is INFERRED from the pair numbering following the
  article's montage rule (sEEG: deep-surface neighbours, adjacent numbers; ECoG grids: row-by-row neighbours, non-adjacent
  numbers); the release itself does not state the electrode type per channel.
- `sub-NN_space-Other_electrodes.tsv`: `data.elec.chanpos` of each bipolar channel as released. The release does not state
  whether chanpos is in native or MNI space (the article localizes in native space and transforms to MNI for group display),
  hence `space-Other`; see coordsystem.json. The per-subject `yeo.mat` (Yeo-7 network per task-responsive electrode) lists
  fewer entries than channels and carries no channel names, so it could not be mapped to channels; it is kept in sourcedata/
  (counts per subject in code/conversion_report.json, key yeo_unmapped). Where its length matches the channel count the
  network is written to the `yeo7_network` column.
- `participants.tsv`: no per-subject age or sex in the release (article Table S1).
- `sourcedata/osf-rx2zd/`: every release file (data incl. `behavior.mat`, `trialinfo.mat`, scripts) unchanged apart from MAT
  text-header dates reduced to month/year; `DEIDENTIFICATION_MANIFEST.tsv` lists sha-256 before and after.

## Privacy
Subjects are numbered S1-S11; files hold trial parameters, responses and reaction times only.
""",
    changes_note="BIDS (derivative) repackaging of the OSF rx2zd release (Johnson et al. 2023) for NEMAR (o5c lane A2).")
json.dump(report, open(os.path.join(OUT, "code", "conversion_report.json"), "w"), indent=1, default=str)
print(json.dumps({"recordings": len(Rr), "subjects": len(subs), "hours": hrs, "rt_all": all(r["roundtrip_ok"] for r in Rr), "channels": [r["n_channels"] for r in Rr], "types": [r["types"] for r in Rr]}, indent=1))
assert all(r["roundtrip_ok"] for r in Rr)
print("CONVERT_DONE")
