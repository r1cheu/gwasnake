def align_to_bed(bed, ids):
    pos = {iid: i for i, iid in enumerate(bed.iid)}
    present = [iid for iid in bed.iid if iid in ids]
    return present, [pos[iid] for iid in present]
