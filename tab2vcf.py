#!/usr/bin/env python
import argparse
import pandas as pd


class Variant(object):
    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.samples = []
        self.sample_info = {}

    def __str__(self):
        return "<Variant @ chr {0}, {1}, ref={2}, alt={3}, samples={4}>".format(
            self.chrom,
            self.pos,
            self.ref,
            self.alt,
            self.samples)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return (self.chrom == other.chrom) and \
                (self.pos == self.pos) and \
                (self.ref == self.ref) and \
                (self.alt == self.alt)

    def add_sample(self, new_sample):
        if new_sample not in self.samples:
            self.samples.append(new_sample)

    def add_info(self, sample, alt_reads, tot_reads):
        if sample in self.sample_info.keys():
            print("WARN: adding the same sample twice to variant record")
        self.sample_info[sample] = (alt_reads, tot_reads)

    def present_in_sample(self, s):
        return s in self.sample_info.keys()

    def get_sample_info(self, s):
        return self.sample_info[s]


def fix_indices(df, colname="SampleID"):
    """
        the samples in the raw data frame
        are in the same column as patient ids,
        this function splits them into two different columns
    """
    def cleaner(value):
        if value[0:2] == "M_":
            parts = value.split("_")
            return parts[-1]
        return value

    split_columns = df[colname].str.split(pat="_", n=1, expand=True)
    split_columns["sample_id"] = split_columns.iloc[:, 1].map(cleaner)
    clean_names = split_columns.drop(1, axis=1)
    clean_names.columns = ["patient_id", "sample_id"]
    return clean_names


def process_patient(name, grp):
    variants = []
    for i, row in grp.iterrows():
        v = Variant(row.Chr, row.Position, row.Ref, row.Alt)
        v.add_sample(row.sample_id)
        v.add_info(row.sample_id, row.AltCount, row.RefCount)
        if v in variants:
            v_ref = variants[variants.index(v)]
            v_ref.add_sample(row.sample_id)
            v_ref.add_info(row.sample_id, row.AltCount, row.RefCount)
        else:
            variants.append(v)
    return variants


def generate_vcf_file(name, variants, all_samples, outdir):
    header = "##fileformat=VCFv4.0\n"
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    all_samples = sorted(all_samples)
    for sample in all_samples:
        header += "\t{0}".format(sample)
    header += "\n"

    for v in variants:
        header += "{0}\t{1}\t.\t{2}\t{3}\t-1\tPASS\t".format(
            v.chrom,
            v.pos,
            v.ref,
            v.alt)
        header += "SOMATIC\tGT:AD:DP\t"
        for sample in all_samples:
            if v.present_in_sample(sample):
                alt_reads, tot_reads = v.get_sample_info(sample)
                header += "0/1:{0}:{1}\t".format(
                    alt_reads,
                    alt_reads + tot_reads)
            else:
                header += "0/0:.:.\t"
        header += "\n"

    output_file = "{0}/{1}.vcf".format(outdir, name)
    with open(output_file, "w") as handle:
        handle.write(header)


def main():
    parser = argparse.ArgumentParser(description='convert .xlsx files to vcf.')
    parser.add_argument('input_file',
                        nargs=1,
                        help="The input file (.xlsx-format)")
    parser.add_argument('vcf_out',
                        nargs=1,
                        default="vcfs/",
                        help="output directory")

    args = parser.parse_args()

    source = pd.read_excel(args.input_file)
    split_columns = fix_indices(source)
    cleaned_data = source.drop("SampleID", axis=1)
    dataset = pd.concat([split_columns, cleaned_data], axis=1)
    dataset.set_index("patient_id")

    grouped = dataset.groupby("patient_id")
    for name, grp in grouped:
        variants = process_patient(name, grp)
        generate_vcf_file(name,
                          variants,
                          list(set(grp.sample_id.tolist())),
                          args.vcf_out)


if __name__ == "__main__":
    main()
